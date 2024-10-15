#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
int main()
{
Eigen::VectorXd v1(3);
Eigen::VectorXd v2(3);
v1 << 1, 2, 3;
v2 << 4, 5, 6;
std::cout << "v1 + v2 = " << v1 + v2 << std::endl;

}
