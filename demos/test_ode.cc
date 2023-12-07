#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;
using namespace std;


class MassSpring : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    // cout << "EvaluateDeriv" << endl;
    // cout << "df = " << endl << df << endl;
    df(0,1) = 1;
    df(1,0) = -1;

    df(0, 0) = 0;
    df(1, 1) = 0;

    // cout << "df = " << endl << df << endl;
  }
};


int main()
{
  double tend = 4*M_PI;
  int steps = 100;
  cout << "tend = " << tend << ", steps = " << steps << endl;
  Vector<> y { 1, 0 };
  cout << "y = " << y << endl;
  auto rhs = make_shared<MassSpring>();
  cout << "rhs = " << rhs << endl;
  

  cout << "SolveODE_IE" << endl;
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });

  cout << "y = " << y << endl;


  Vector<> y_CR { 1, 0 };
  cout << "SolveODE_CrankNicolson" << endl;
  SolveODE_CrankNicolson(tend, steps, y_CR, rhs,
                         [](double t, VectorView<double> y_CR) { cout << t << "," << y_CR(0) << "," << y_CR(1) << endl; });

  cout << "y = " << y << endl;

  // Exiting this method throws a double free error. I have no idea why, free() is never used.
  return 0;
}
