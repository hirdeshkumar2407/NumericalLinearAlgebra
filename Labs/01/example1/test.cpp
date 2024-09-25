#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
int main()
{
  // Create a 3x3 matrix with random values
  MatrixXd m = MatrixXd::Random(3,3);

  // Add a constant matrix (1.0) to the random matrix and multiply by 10
  m = (m + MatrixXd::Constant(3,3,1.0)) * 10;
  std::cout << "m =" << std::endl << m << std::endl;
 
     // Create a vector with 3 elements
  VectorXd v(3);

   // Output the product of matrix m and vector v
  v << 1, 0, 0;
  std::cout << "m * v =" << std::endl << m * v << std::endl;
}