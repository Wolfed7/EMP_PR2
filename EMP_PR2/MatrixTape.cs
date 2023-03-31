namespace EMP_PR2;
public class MatrixTape
{
   private double[] Diag;
   private double[][] Lower;
   private double[][] Upper;
   public int TapeWidth { get; init; }
   public int DiagLength { get; init; }

   public MatrixTape(int tapeWidth, int diagLength)
   {
      TapeWidth = tapeWidth;
      DiagLength = diagLength;

      Diag = new double[DiagLength];
      Upper = new double[diagLength].Select(row => new double[(tapeWidth - 1) / 2]).ToArray();
      Lower = new double[diagLength].Select(row => new double[(tapeWidth - 1) / 2]).ToArray();
   }
}

