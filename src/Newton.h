#ifndef Newton_h
#define Newton_h

#include "nonlinfunc.h"

namespace ASC_ode
{

  void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<> res(func->DimF());
    Matrix<> fprime(func->DimF(), func->DimX());

    for (int i = 0; i < maxsteps; i++)
      {
        cout << "Newton step " << i << endl;
        func->Evaluate(x, res);
        cout << "res = " << endl << res << endl;
        // cout << "|res| = " << L2Norm(res) << endl;
        fprime *= 0;
        cout << "fprime = " << endl << fprime << endl;
        func->EvaluateDeriv(x, fprime);
        cout << "fprime = " << endl << fprime << endl;
        //cout << "fprime^-1 = " << endl << Matrix(InverseLapack(fprime)) << endl;
        //cout << "fprime^-1*fprime = " << endl << Matrix(Matrix(fprime*InverseLapack(fprime))) << endl;
        fprime = Matrix(InverseLapack(fprime));
        cout << "fprime = " << endl << fprime << endl;
        x -= fprime*res;
        cout << "frprime*res = " << endl << fprime*res << endl;
        cout << "x = " << endl << x << endl;

        double err= res.Norm2();
        if (callback)
          callback(i, err, x);

        cout << "err = " << err << endl;

        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
