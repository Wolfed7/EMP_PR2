namespace EMP_PR2;

public class Matrix
{
   private double[,] _container;
   public int Size { get; }

   public double this[int i, int j]
   {
      get => _container[i, j];
      set => _container[i, j] = value;
   }

   public Matrix(int size)
   {
      _container = new double[size, size];
      Size = size;
   }

   public void Clear()
       => Array.Clear(_container, 0, _container.Length);

   public void Copy(Matrix destination)
   {
      for (int i = 0; i < destination.Size; i++)
         for (int j = 0; j < destination.Size; j++)
            destination[i, j] = _container[i, j];
   }

   public static Matrix operator +(Matrix fstMatrix, Matrix sndMatrix)
   {
      Matrix resultMatrix = new(fstMatrix.Size);

      for (int i = 0; i < resultMatrix.Size; i++)
         for (int j = 0; j < resultMatrix.Size; j++)
            resultMatrix[i, j] = fstMatrix[i, j] + sndMatrix[i, j];

      return resultMatrix;
   }

   public static Matrix operator *(double coef, Matrix Mat)
   {
      Matrix resultMatrix = new(Mat.Size);

      for (int i = 0; i < resultMatrix.Size; i++)
         for (int j = 0; j < resultMatrix.Size; j++)
            resultMatrix[i, j] = coef * Mat[i, j];

      return resultMatrix;
   }
}