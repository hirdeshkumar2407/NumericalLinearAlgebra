#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
 
int main()
{
   MatrixXd A = MatrixXd::Random(100,100);
   MatrixXd b = MatrixXd::Random(100,50);
   // HERE WE ARE SOLVING 50 SYSTEMS OF EQUATION
   MatrixXd x = A.fullPivLu().solve(b);
   double relative_residual = (A*x - b).norm() / b.norm(); // norm() is L2 norm EUCLDEAN NORM
   std::cout << "The relative residual is:\n" << relative_residual << std::endl;
   //since matrix is random, it is not guranted you get good solution
}
