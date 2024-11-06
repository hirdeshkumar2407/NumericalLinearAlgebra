#include <cstdlib> // System includes
#include <iostream>
#include <Eigen/SparseCore> // Eigen includes
#include <Eigen/IterativeLinearSolvers> // Eigen includes
#include <unsupported/Eigen/SparseExtra>

using std::endl;
using std::cout;
#include "cg.hpp" // Local includes

/*We test the conjugate gradient method on a tri-diagonal linear system assuming that the exact solution as all the coefficients equal */

int main(int argc, char** argv)
{

    using namespace LinearAlgebra; // Using namespace for convenience

    using SpMat = Eigen::SparseMatrix<double>; // Alias for SparseMatrix
    using SpVec = Eigen::VectorXd; // Alias for VectorXd

    int n = 1000; // Size of the matrix
    SpMat A(n, n); // Define matrix
    A.reserve(2998); // Reserve space for non-zero entries
    for (int i = 0; i <n; i++){
        A.coeffRef(i,i ) = 2.0 * (i + 1); // Diagonal entries
        if (i > 0) A.coeffRef(i, i - 1) = -i; // Lower diagonal entries
        if (i < n - 1) A.coeffRef(i, i + 1) = -(i + 1); // Upper diagonal entries   
    }
    /*we set the parameters for the conjugate gradient routine (desired tolerance, initial guess, and maximum number of iterations allowed for convergence). As precoditioner we take the simple diagonal preconditioner provided by Eigen (also called the Jacobi preconditioner).*/

    double tol = 1.e-10; // Convergence tolerance

    int result, maxit = 1000; // Maximum iterations

   std::string matrixFileOut = "./AtestCG.mtx";
    Eigen::saveMarket(A, matrixFileOut);

    std::cout << "Matrix A has been saved to " << matrixFileOut << std::endl;
    std::cout << "Matrix size: " << A.rows() << "x" << A.cols() << std::endl;
    std::cout << "Non-zero entries: " << A.nonZeros() << std::endl;

    SpMat B = SpMat(A.transpose()) - A; // Check symmetry
    std::cout << "Norm of A - A.t: " << B.norm() << std::endl;

    // Create Rhs b
    SpVec e = SpVec::Ones(A.rows());
    SpVec b = A * e;
    SpVec x(A.rows());

    // Eigen Jacobi

    
    Eigen::DiagonalPreconditioner<double> D(A); // Create diagonal preconditioner

    // First with eigen CG

    Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
    cg.setMaxIterations(maxit);
    cg.setTolerance(tol);
    cg.compute(A);
    x = cg.solve(b);
    std::cout << " Eigen native CG" << std::endl;
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    std::cout << "effective error: " << (x - e).norm() << std::endl;

     // Now with hand-made CG
  x=0*x;
  result = CG(A, x, b, D, maxit, tol);        // Solve system

  std::cout <<" hand-made CG "<< endl;
  cout << "CG flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;

  return result;

}
