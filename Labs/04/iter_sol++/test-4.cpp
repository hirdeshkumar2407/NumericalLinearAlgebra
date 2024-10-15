
#include "cgs.hpp"
#include "bcgstab.hpp"
#include "gmres.hpp"

int main(int argc, char** argv)
{
  using namespace LinearAlgebra;
  // Some useful alias  ...

  int n = 1000;
  double gam = -0.5;
  SpMat A(n,n);                      // define matrix
  A.reserve(2997);
  for (int i=0; i<n; i++) {
      A.coeffRef(i, i) = 2.0;
      if(i>1) A.coeffRef(i, i-2) = gam;
      if(i<n-1) A.coeffRef(i, i+1) = 1.0;
  }

  double tol = 1.e-8;                // Convergence tolerance
  int result, maxit = 1000;           // Maximum iterations
  int restart = maxit;                  // Restart for gmres

  std::cout<<"Matrix size:"<<A.rows()<<"X"<<A.cols()<<std::endl;
  std::cout<<"Non zero entries:"<<A.nonZeros()<<std::endl;
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());
  Eigen::LeastSquareDiagonalPreconditioner<double> SD(A);

  // Solve with CGS method
  x=0*x;
  result = CGS(A, x, b, SD, maxit, tol);
  cout << "CGS: #iterations: " << result << endl;

  // Solve with BiCGSTAB method
  ...

  // Solve with GMRES method
  ...

  return 0;
}