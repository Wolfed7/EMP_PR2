namespace EMP_PR2;
public class MatrixTape
{
   public double[] Diag;
   public double[][] Lower;
   public double[][] Upper;
   public int TapeWidth { get; init; }
   public int DiagLength { get; init; }
   public int TriangleWidth { get; init; }

   public MatrixTape(int tapeWidth, int diagLength)
   {
      TapeWidth = tapeWidth;
      DiagLength = diagLength;
      TriangleWidth = (TapeWidth - 1) / 2;

      Diag = new double[DiagLength];
      Upper = new double[DiagLength].Select(row => new double[TriangleWidth]).ToArray();
      Lower = new double[DiagLength].Select(row => new double[TriangleWidth]).ToArray();
   }

   public void Clear()
   {
      Array.Clear(Diag);

      foreach (var item in Lower)
      {
         Array.Clear(item);
      };

      foreach (var item in Upper)
      {
         Array.Clear(item);
      };
   }
}

