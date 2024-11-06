#include <cstdlib>                      // System includes
#include <iostream>
#include <Eigen/SparseCore>             // Eigen includes
#include <Eigen/IterativeLinearSolvers> // Eigen includes

using std::endl;
using std::cout;
#include "cg.hpp" // Local includes

/* We test the Gradient method (Steepest Descent) on a tri-diagonal linear system. 
   The exact solution is assumed to be all coefficients equal. */

namespace LinearAlgebra
{
template <class Matrix, class Vector>
int GradientMethod(const Matrix &A, Vector &x, const Vector &b, int &max_iter, typename Vector::Scalar &tol)
{
    using Real = typename Matrix::Scalar;
    Real resid;

    // Step 1: Calculate initial residual
    Vector r = b - A * x;
    Vector p = r; // Steepest descent direction is just the residual initially
    Real normb = b.norm();
    
    if (normb == 0.0) normb = 1;

    // Check if initial guess is good enough
    if ((resid = r.norm() / normb) <= tol)
    {
        tol = resid;
        max_iter = 0;
        return 0; // Converged immediately
    }

    // Step 2: Main loop of the Gradient Method
    for (int i = 1; i <= max_iter; ++i)
    {
        // Compute alpha (step size)
        Real alpha = r.dot(r) / (p.dot(A * p));
        
        // Update solution x
        x += alpha * p;

        // Update residual
        r = b - A * x;

        // Check for convergence
        if ((resid = r.norm() / normb) <= tol)
        {
            tol = resid;
            max_iter = i;
            return 0; // Converged
        }
    }

    // If not converged within max iterations
    tol = resid;
    return 1; // Failed to converge
}
} // namespace LinearAlgebra

int main(int argc, char **argv)
{
    using namespace LinearAlgebra;

    using SpMat = Eigen::SparseMatrix<double>;
    using SpVec = Eigen::VectorXd;

    int n = 1000; // Size of the matrix
    SpMat A(n, n); // Define matrix
    A.reserve(2998); // Reserve space for non-zero entries
    for (int i = 0; i < n; i++)
    {
        A.coeffRef(i, i) = 2.0 * (i + 1); // Diagonal entries
        if (i > 0)
            A.coeffRef(i, i - 1) = -i; // Lower diagonal entries
        if (i < n - 1)
            A.coeffRef(i, i + 1) = -(i + 1); // Upper diagonal entries
    }

    // Set the parameters for the gradient method
    double tol = 1.e-10; // Convergence tolerance
    int result, maxit = 1000; // Maximum iterations

    std::cout << "Matrix size: " << A.rows() << "x" << A.cols() << std::endl;
    std::cout << "Non-zero entries: " << A.nonZeros() << std::endl;

    SpMat B = SpMat(A.transpose()) - A; // Check symmetry
    std::cout << "Norm of A - A.t: " << B.norm() << std::endl;

    // Create Rhs b
    SpVec e = SpVec::Ones(A.rows());
    SpVec b = A * e;
    SpVec x(A.rows()); // Initial guess
    x.setZero(); // Initial guess is zero

    // Now with Gradient Method (Steepest Descent)
    result = GradientMethod(A, x, b, maxit, tol); // Solve system

    std::cout << "Gradient Method (Steepest Descent)" << std::endl;
    std::cout << "Result = " << result << std::endl;
    std::cout << "Iterations performed: " << maxit << std::endl;
    std::cout << "Effective error: " << (x - e).norm() << std::endl;

    return result;
}
