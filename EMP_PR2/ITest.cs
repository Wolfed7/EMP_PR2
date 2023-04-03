namespace EMP_PR2;

public interface ITest
{
   public Func<double, double> F(Func<double, double> f); // Зависит от U
   public double Lambda(double x, double t);
   public double Sigma(double x, double t);

   // Параметр добавлен в целях учебной задачи для исследований,
   // программа строится на нём и не имеет способности
   // нормального задания реальной задачи.
   public double U(double x, double t);

   public double DerivativeF(double u);
}

// u = x + t
public class Test1 : ITest
{
   public Func<double, double> F(Func<double, double> f)
      => (u) => f(u);

   public double Lambda(double x, double t)
      => 1;

   public double Sigma(double x, double t)
      => x + t;

   public double U(double x, double t)
      => x + t;

   public double DerivativeF(double u)
      => 1;
}

// u = x * t
public class Test2 : ITest
{
   public double U(double x, double t)
      => x * t;

   public Func<double, double> F(Func<double, double> f)
      => (u) => f(u) * f(u);

   public double Lambda(double x, double t)
      => x;

   public double Sigma(double x, double t)
      => x * t * t + t / x;

   public double DerivativeF(double u)
   => 2 * u;
}

// x^2 * t
public class Test3 : ITest
{
   public Func<double, double> F(Func<double, double> f)
      => (u) => f(u);

   public double Lambda(double x, double t)
      => 1;

   public double Sigma(double x, double t)
      => t + 2 * t / x / x;

   public double U(double x, double t)
      => x * x + t;

   public double DerivativeF(double u)
      => 1;
}

// u = t * cos(x)
public class Test4 : ITest
{
   public Func<double, double> F(Func<double, double> f)
      => (u) => Math.Exp(f(u));

   public double Lambda(double x, double t)
      => 1;

   public double Sigma(double x, double t)
      => Math.Exp(t * Math.Cos(x)) / Math.Cos(x) - t;

   public double U(double x, double t)
      => t * Math.Cos(x);

   public double DerivativeF(double u)
      => Math.Exp(u);
}

// u = t * ch(x)
public class Test5 : ITest
{
   public Func<double, double> F(Func<double, double> f)
      => (u) => Math.Exp(f(u));

   public double Lambda(double x, double t)
      => 1;

   public double Sigma(double x, double t)
      => Math.Exp(t * Math.Cosh(x)) / Math.Cosh(x) + t;

   public double U(double x, double t)
      => t * Math.Cosh(x);

   public double DerivativeF(double u)
      => Math.Exp(u);
}

