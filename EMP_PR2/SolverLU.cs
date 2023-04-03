namespace EMP_PR2;

public class SolverLU
{
   private double[] RightPart;
   private MatrixTape Matrix;
   private double[] Result;

   public SolverLU(MatrixTape matrix, double[] rightPart)
   {
      RightPart = rightPart;
      Matrix = matrix;
      Result = new double[rightPart.Length];
   }

   public double[] Compute()
   {
      Result = new double[RightPart.Length];
      Array.Copy(RightPart, Result, RightPart.Length);
      LUDecompose();
      Forward();
      Backward();
      return Result;
   }

   // LU разложение для матрицы в ленточном формате (с отдельной диагональю)
   private void LUDecompose()
   {
      for (int i = 0; i < Matrix.DiagLength; i++)
      {

         for (int j = 0; j < Matrix.TriangleWidth; j++)
         {
            int j0 = i + j - Matrix.TriangleWidth;
            if (j0 < 0)
               continue;

            double sumL = 0.0;
            double sumU = 0.0;

            for (int ik = 0; ik < j; ik++)
            {
               int jk = Matrix.TriangleWidth - j + ik;
               sumL += Matrix.Lower[i][ik] * Matrix.Upper[j0][jk];
               sumU += Matrix.Upper[i][ik] * Matrix.Lower[j0][jk];
            }

            Matrix.Lower[i][j] -= sumL;
            Matrix.Upper[i][j] -= sumU;
            Matrix.Upper[i][j] /= Matrix.Diag[j0];
         }

         double sumD = 0.0;
         for (int j = 0; j < Matrix.TriangleWidth; j++)
            sumD += Matrix.Lower[i][j] * Matrix.Upper[i][j];

         Matrix.Diag[i] -= sumD;
      }
   }

   // Прямой обход по слау.
   private void Forward()
   {
      for (int i = 0; i < Matrix.DiagLength; i++)
      {
         for (int j = 0; j < Matrix.TriangleWidth; j++)
         {
            int j0 = i - Matrix.TriangleWidth + j;
            if (j0 < 0)
               continue;

            Result[i] -= Matrix.Lower[i][j] * Result[j0];
         }

         Result[i] /= Matrix.Diag[i];
      }
   }

   // Обратный обход по слау.
   private void Backward()
   {
      for (int i = Matrix.DiagLength - 1; i >= 0; i--)
      {
         for (int j = Matrix.TriangleWidth - 1; j >= 0; j--)
         {
            int i0 = i - Matrix.TriangleWidth + j;
            if (i0 < 0)
               continue;

            Result[i0] -= Matrix.Upper[i][j] * Result[i];
         }
      }
   }
}
