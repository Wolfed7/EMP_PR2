namespace EMP_PR2;

public class Mesh
{
   public double[] Nodes { get; private set; }
   public int Length => Nodes.Length;
   public double StartPoint { get; init; } // Левая граница отрезка.
   public double EndPoint { get; init; } // Правая граница отрезка.
   private int SectionsCount; // Количество разбиений сетки.
   private double DischargeRatio; // Коэффициент разрядки сетки.

   public Mesh(string spaceMeshPath)
   {
      try
      {
         using (StreamReader sr = new(spaceMeshPath))
         {
            StartPoint = double.Parse(sr.ReadLine());
            EndPoint = double.Parse(sr.ReadLine());

            SectionsCount = int.Parse(sr.ReadLine());
            DischargeRatio = double.Parse(sr.ReadLine());

            Nodes = new double[SectionsCount + 1];
         }
      }
      catch (Exception)
      {

         throw;
      }

      BuildMesh();
   }

   //public SpaceMesh(double startPoint, double endPoint, int sectionsCount, int dischargeRatio)
   //{
   //   StartPoint = startPoint;
   //   EndPoint = endPoint;
   //   SectionsCount = sectionsCount;
   //   DischargeRatio = dischargeRatio;
   //   Nodes = new double[sectionsCount + 1];
   //}

   public double this[int i]
   {
      get => Nodes[i];
      set => Nodes[i] = value;
   }

   public void BuildMesh()
   {
      double sum = 0;
      double dischargeRatioPow = 1;
      for (int i = 0; i < SectionsCount; i++)
      {
         dischargeRatioPow *= DischargeRatio;
         sum += dischargeRatioPow;
      }

      double step = (EndPoint - StartPoint) / sum;
      Nodes[0] = StartPoint;
      for (int i = 1; i < SectionsCount; i++)
      {
         Nodes[i] = Nodes[i - 1] + step;
         step *= DischargeRatio;
      }
      Nodes[^1] = EndPoint;
   }
}

