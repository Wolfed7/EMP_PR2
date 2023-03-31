using System.Collections.Immutable;

namespace EMP_PR2;

public abstract class SpaceMesh
{
   public double StartPoint { get; init; } // Левая граница отрезка.
   public double EndPoint { get; init; } // Правая граница отрезка.
   public int SectionsCount { get; init; } // Количество разбиений сетки.

   public ImmutableArray<(int, double, double)> SubAreas;
}

