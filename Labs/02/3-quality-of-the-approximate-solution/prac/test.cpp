#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

int main(){

    MatrixXd A = MatrixXd::Random(100,100);
    MatrixXd b = MatrixXd::Random(100,50);
    /*If you know more about the properties of your matrix, you can use the above table to select the best method. For example, a good choice for solving linear systems with a non-symmetric matrix of full rank is `PartialPivLU`. If you know that your matrix is also symmetric and positive definite, the above table says that a very good choice is the LLT or LDLT decomposition. */
    MatrixXd x = A.fullPivLu().solve(b);

    double relative_resdiual = (A*x - b).norm() / b.norm();
    std::cout << "The relative residual is: " << relative_resdiual << std::endl;
}