namespace EMP_PR2;

public enum NonlinearMethod
{
   SIMPLE_ITER,
   NEWTON
}

public class FEM
{
   private double norm;
   private List<(double, double, double, int)> relaxTest;

   private Mesh _space;
   private Mesh _time;
   private ITest _test;

   // Локальные базисные функции (линейные).
   private delegate double Basis(double x, int ielem);
   private Basis[] _basis;
   private Basis[] _dBasis;

   // Два временных слоя, содержат вектор q,
   // матрица 2 * q
   private double[][] _layers;

   private NonlinearMethod _method;

   // Точность для выхода.
   public double Eps { get; init; }

   // Итерации до аварийного выхода.
   public int MaxIters { get; init; }

   private MatrixTape _globalMatrix;
   private double[] _globalRightPart;

   private double[] _localRightPart;
   private Matrix _localStiffness;
   private Matrix _localMass;

   private SolverLU _SLAE;

   public FEM(Mesh space, Mesh time, ITest test, NonlinearMethod method, double eps, int maxiters)
   {
      _space = space;
      _time = time;
      _test = test;
      _method = method;
      _layers = new double[2][].Select(q => new double[_space.Length]).ToArray();

      _localStiffness = new(2);
      _localMass = new(2);
      _localRightPart = new double[2];

      _globalRightPart = new double[_space.Length];
      _globalMatrix = new(3, _space.Length);
      _basis = new Basis[] { Psi1, Psi2 };
      _dBasis = new Basis[] { dPsi1, dPsi2 };

      Eps = eps;
      MaxIters = maxiters;


      relaxTest = new List<(double, double, double, int)>(0);
   }

   public void Update(ITest test)
      => _test = test;

   public void Compute()
   {
      // Цикл для исследования на лучший параметр релаксации.
      // БЕЗ РЕЛАКСАЦИИ ПОСТАВИТЬ 100 101
      for (int irelax = 70; irelax < 156; irelax++)
      {
         double timeDifference;
         double residual = 0.0;
         double relax = irelax / 100.0;
         int iter;
         int globalIters = 0;
         double[] qPrev = new double[_globalMatrix.DiagLength];

         AccountStartCondition();
         Array.Copy(_layers[0], _layers[1], _layers[0].Length);

         for (int itime = 1; itime < _time.Length; itime++)
         {
            timeDifference = _time[itime] - _time[itime - 1];
            AssemblySlae(itime, timeDifference);
            AccountFirstBoundaryConditions(itime);

            if (_method == NonlinearMethod.NEWTON)
            {
               Linearization();
               AccountFirstBoundaryConditions(itime);
            }

            for (iter = 1; iter < MaxIters; iter++)
            {
               _SLAE = new(_globalMatrix, _globalRightPart);
               _layers[1] = AddVectors(MultVector(relax, _SLAE.Compute()), MultVector(1 - relax, qPrev));

               AssemblySlae(itime, timeDifference);
               AccountFirstBoundaryConditions(itime);

               residual = NormFEM(DiffVectors(MultMatrixOnVec(_globalMatrix, _layers[1]), _globalRightPart)) / NormFEM(_globalRightPart);

               if (_method == NonlinearMethod.NEWTON)
               {
                  Linearization();
                  AccountFirstBoundaryConditions(itime);
               }

               if (NormFEM(DiffVectors(_layers[1], qPrev)) / NormFEM(_layers[1]) < Eps)
                  break;

               Array.Copy(_layers[1], qPrev, _layers[1].Length);

               if (residual < Eps)
                  break;
            }

            globalIters += iter;
            Array.Copy(_layers[1], _layers[0], _layers[1].Length);
         }
         

         // Точное решение для сравнения с полученным.
         double[] exactSolution = new double[_space.Length];
         for (int i = 0; i < exactSolution.Length; i++)
         {
            exactSolution[i] = _test.U(_space[i], _time[^1]);
         }


         //Console.WriteLine
         //(
         //   "Сетка".PadRight(20) +
         //   "Точное решение".PadRight(20) +
         //   "Численное решение".PadRight(20)
         //);

         //for (int i = 0; i < exactSolution.Length; i++)
         //{
         //   Console.Write    ("{0:e3}".PadRight(20), _space[i]);
         //   Console.Write    ("{0:e15}".PadRight(20), exactSolution[i]);
         //   Console.WriteLine("{0:e15}".PadRight(20), _layers[1][i]);
         //}

         norm = Norm(DiffVectors(exactSolution, _layers[1])) / Norm(exactSolution);
         //Console.WriteLine($"Коэффициент релаксации = {relax}");
         //Console.WriteLine($"Относительная погрешность = {norm:0.00E+0}");
         //Console.WriteLine($"Невязка = {residual:0.00E+0}");
         //Console.WriteLine($"Количество итераций = {globalIters}");

         relaxTest.Add( (relax, norm, residual, globalIters) );
      }

      for (int i = 0; i < relaxTest.Count; i++)
      {
         Console.Write("{0:f2}".PadRight(20), relaxTest[i].Item1);
         Console.Write("{0:e2}".PadRight(20), relaxTest[i].Item2);
         Console.Write("{0:e2}".PadRight(20), relaxTest[i].Item3);
         Console.WriteLine("{0:f0}".PadRight(20), relaxTest[i].Item4);
      }
   }

   private double NormFEM(double[] vector)
   {
      double result = 0;

      // Избегаем краевых в точках 0 и vector.Length - 1.
      for (int i = 1; i < vector.Length - 1; i++)
         result += vector[i] * vector[i];

      return Math.Sqrt(result);
   }

   // Учёт начальных условий.
   // Не работает при реальных задачах, учебный костыль.
   private void AccountStartCondition()
   {
      for (int i = 0; i < _layers[0].Length; i++)
         _layers[0][i] = _test.U(_space.Nodes[i], _time.Nodes[0]);
   }

   // Учёт первых краевых условий.
   // Они тут в принципе одни, лаба не про это.
   private void AccountFirstBoundaryConditions(int itime)
   {
      _globalMatrix.Diag[0] = 1;
      _globalMatrix.Diag[^1] = 1;

      _globalMatrix.Upper[1][0] = 0;
      _globalMatrix.Lower[^1][0] = 0;

      _globalRightPart[0]  = _test.U(_space[0] , _time[itime]);
      _globalRightPart[^1] = _test.U(_space[^1], _time[itime]);
   }

   private void AssemblySlae(int itime, double timeDiff)
   {
      Array.Clear(_globalRightPart);
      _globalMatrix.Clear();

      for (int ielem = 0; ielem < _space.Length - 1; ielem++)
      {

         AssemblyLocalMatrixes(ielem, itime, timeDiff);
         AssemblyLocalRightPart(ielem);

         // Получим локальную матрицу A
         // (Сумма локальных матрицы жёсткости и матрицы массы).
         _localStiffness += _localMass;

         for (int i = 0; i < _localStiffness.Size; i++)
            for (int j = 0; j < _localStiffness.Size; j++)
               AddElemToGlobalMatrix(ielem + i, ielem + j, _localStiffness[i, j]);

         AddElemToGlobalRightPart(ielem);

         _localStiffness.Clear();
         _localMass.Clear();
         Array.Clear(_localRightPart);
      }
   }

   private void AssemblyLocalMatrixes(int ielem, int itime, double timeDiff)
   {
      for (int i = 0; i < _localStiffness.Size; i++)
         for (int j = 0; j < _localStiffness.Size; j++)
            _localStiffness[i, j] = GaussEdge(_test.Lambda, _dBasis[i], _dBasis[j], ielem, itime);

      for (int i = 0; i < _localMass.Size; i++)
         for (int j = 0; j < _localMass.Size; j++)
            _localMass[i, j] = GaussEdge(_test.Sigma, _basis[i], _basis[j], ielem, itime) / timeDiff;
   }

   private void AssemblyLocalRightPart(int ielem)
   {
      // Аппроксимация функции на элементе. uh = q * Psi
      double Uh(double x)
         => _layers[1][ielem] * _basis[0](x, ielem)
         + _layers[1][ielem + 1] * _basis[1](x, ielem);

      for (int i = 0; i < _localRightPart.Length; i++)
      {
         // Здесь сборка вектора b1 численным интегрированием.
         _localRightPart[i] = GaussEdge(_test.F(Uh), _basis[i], ielem);

         // Здесь досборка вектора d = b1 + 1/delta t * M * q0 (q0 - вектор весов на предыдущем слое)
         // (1/delta t уже в локальной матрице масс).
         for (int j = 0; j < _localMass.Size; j++)
            _localRightPart[i] += _layers[0][ielem + j] * _localMass[i, j];
      }
   }

   private void AddElemToGlobalMatrix(int i, int j, double value)
   {
      if (i == j)
         _globalMatrix.Diag[i] += value;
      else if (i > j)
         _globalMatrix.Lower[i][0] += value;
      else
         _globalMatrix.Upper[j][0] += value;
   }

   private void AddElemToGlobalRightPart(int ielem)
   {
      for (int i = 0; i < _localRightPart.Length; i++)
         _globalRightPart[ielem + i] += _localRightPart[i];
   }

   private double[] MultMatrixOnVec(MatrixTape matrix, double[] vector)
   {
      double[] result = new double[vector.Length];

      for (int i = 0; i < vector.Length; i++)
      {
         result[i] += matrix.Diag[i] * vector[i];
         for (int j = 0; j < matrix.TriangleWidth; j++)
         {
            int j0 = i + j - matrix.TriangleWidth;
            if (j0 < 0)
               continue;

            result[i] += matrix.Lower[i][j] * vector[j0];
            result[j0] += matrix.Upper[i][j] * vector[i];
         }
      }

      return result;
   }

   private double[] MultVector(double value, double[] vector)
   {
      double[] newvec = new double[vector.Length];
      for (int i = 0; i < vector.Length; i++)
         newvec[i] = value * vector[i];
      return newvec;
   }

   private double[] AddVectors(double[] vector1, double[] vector2)
   {
      double[] sum = new double[vector1.Length];
      for (int i = 0; i < vector1.Length; i++)
         sum[i] = vector1[i] + vector2[i];
      return sum;
   }

   private double[] DiffVectors(double[] vector1, double[] vector2)
   {
      double[] sum = new double[vector1.Length];
      for (int i = 0; i < vector1.Length; i++)
         sum[i] = vector1[i] - vector2[i];
      return sum;
   }

   public double Psi1(double x, int ielem)
   => (_space[ielem + 1] - x) / (_space[ielem + 1] - _space[ielem]);
   public double Psi2(double x, int ielem)
      => (x - _space[ielem]) / (_space[ielem + 1] - _space[ielem]);
   public double dPsi1(double x, int ielem)
      => -1 / (_space[ielem + 1] - _space[ielem]);
   public double dPsi2(double x, int ielem)
      => 1 / (_space[ielem + 1] - _space[ielem]);

   // Пятый Гаусс избыточен (я бы даже сказал вообще лишнее, матрички 2x2 аналитические можно)
   private double GaussEdge(Func<double, double, double> function, Basis psiI, Basis psiJ, int ielem, int itime)
   {
      double fstPoint = _space[ielem];
      double sndPoint = _space[ielem + 1];

      double[] p = 
      { 
         0.0,
         1.0 / 3.0 * Math.Sqrt(5.0 - 2.0 * Math.Sqrt(10.0 / 7.0)),
        -1.0 / 3.0 * Math.Sqrt(5.0 - 2.0 * Math.Sqrt(10.0 / 7.0)),
         1.0 / 3.0 * Math.Sqrt(5.0 + 2.0 * Math.Sqrt(10.0 / 7.0)),
        -1.0 / 3.0 * Math.Sqrt(5.0 + 2.0 * Math.Sqrt(10.0 / 7.0))
      };

      double[] w = 
      { 
         128.0 / 225.0,
         (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
         (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
         (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0,
         (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0 
      };

      double res = 0;

      double lengthEdge = Math.Abs(sndPoint - fstPoint);

      for (int i = 0; i < w.Length; i++)
      {
         double point = (sndPoint - fstPoint) * (1 + p[i]) / 2 + fstPoint;
         res += function(point, _time[itime]) * psiI(point, ielem) * psiJ(point, ielem) * w[i];
      }

      return lengthEdge * res / 2;
   }

   // Этот для сборки вектора правой части, уже более оправдан, но можно интерполянтом.
   private double GaussEdge(Func<double, double> function, Basis psiI, int ielem)
   {
      double fstPoint = _space[ielem];
      double sndPoint = _space[ielem + 1];

      double[] p = 
      { 
         0.0,
         1.0 / 3.0 * Math.Sqrt(5.0 - 2.0 * Math.Sqrt(10.0 / 7.0)),
        -1.0 / 3.0 * Math.Sqrt(5.0 - 2.0 * Math.Sqrt(10.0 / 7.0)),
         1.0 / 3.0 * Math.Sqrt(5.0 + 2.0 * Math.Sqrt(10.0 / 7.0)),
        -1.0 / 3.0 * Math.Sqrt(5.0 + 2.0 * Math.Sqrt(10.0 / 7.0))
      };

      double[] w = 
      { 
         128.0 / 225.0,
         (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
         (322.0 + 13.0 * Math.Sqrt(70.0)) / 900.0,
         (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0,
         (322.0 - 13.0 * Math.Sqrt(70.0)) / 900.0
      };

      double res = 0;

      double lengthEdge = Math.Abs(sndPoint - fstPoint);

      for (int i = 0; i < w.Length; i++)
      {
         double point = (sndPoint - fstPoint) * (1 + p[i]) / 2 + fstPoint;
         res += function(point) * psiI(point, ielem) * w[i];
      }

      return lengthEdge * res / 2;
   }

   private double Norm(double[] vec)
   {
      double result = 0;

      for (int i = 0; i < vec.Length; i++)
         result += vec[i] * vec[i];

      return Math.Sqrt(result);
   }

   private void Linearization()
   {
      double length;
      for (int ielem = 0; ielem < _space.Length - 1; ielem++)
      {
         length = _space[ielem + 1] - _space[ielem];
         _localStiffness[0, 0] = -_test.DerivativeF(_layers[1][ielem]) * length / 3;
         _localStiffness[0, 1] = -_test.DerivativeF(_layers[1][ielem + 1]) * length / 6;
         _localStiffness[1, 0] = -_test.DerivativeF(_layers[1][ielem]) * length / 6;
         _localStiffness[1, 1] = -_test.DerivativeF(_layers[1][ielem + 1]) * length / 3;

         _localRightPart[0] = -_test.DerivativeF(_layers[1][ielem]) * length / 3 * _layers[1][ielem] - _test.DerivativeF(_layers[1][ielem + 1]) * length / 6 * _layers[1][ielem + 1];
         _localRightPart[1] = -_test.DerivativeF(_layers[1][ielem]) * length / 6 * _layers[1][ielem] - _test.DerivativeF(_layers[1][ielem + 1]) * length / 3 * _layers[1][ielem + 1];

         for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
               AddElemToGlobalMatrix(ielem + i, ielem + j, _localStiffness[i, j]);

         AddElemToGlobalRightPart(ielem);

         _localStiffness.Clear();
         _localMass.Clear();
         Array.Clear(_localRightPart);
      }
   }
}
