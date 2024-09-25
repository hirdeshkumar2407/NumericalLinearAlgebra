#include <iostream>
#include <Eigen/Dense>

 
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
int main(){
MatrixXd m(3,3);

m << 1,2,3,
     4,5,6,
     7,8,9;

std::cout << "m matrix" << std::endl << m << std::endl;
std::cout << "m size ="  << m.size() << std::endl;
std::cout << "m rows ="  << m.rows() << std::endl;
std::cout << "m cols ="  << m.cols() << std::endl;
std::cout << "m transpose ="  << std::endl << m.transpose() << std::endl;
std::cout << "m diagonal ="  << std::endl  << m.diagonal() << std::endl;
std::cout << "m adjoint ="  << std::endl  << m.adjoint() << std::endl;
std::cout << "m transpose().colwise().reverse()  ="  << std::endl  << m.transpose().colwise().reverse()   << std::endl;
std::cout << "m rowwise().reverse() ="  << std::endl  << m.rowwise().reverse() << std::endl;
std::cout << "m colwise().reverse()  ="  << std::endl  << m.colwise().reverse()  << std::endl;
std::cout << "m .replicate(i,j)"  << std::endl  << m.replicate(2,2) << std::endl;
printf("\n\n\n\n");

MatrixXd m2= MatrixXd::Identity(3,3);
std::cout << "m2 matrix" << std::endl << m2 << std::endl;
std::cout << "m2 size ="  << m2.size() << std::endl;
std::cout << "m2 rows ="  << m2.rows() << std::endl;
std::cout << "m2 cols ="  << m2.cols() << std::endl;
printf("\n\n\n\n");


MatrixXd m3= MatrixXd::Identity(4,4);
std::cout << "before m4 matrix" << std::endl << m3 << std::endl;
m3.setIdentity(2,2);
std::cout << "after m4 matrix" << std::endl << m3 << std::endl;
printf("\n\n\n\n");

MatrixXd m4(2,3);
m4 << 1,2,3,
     4,5,6;
std::cout << "before resize m5 matrix" << std::endl << m4 << std::endl; 
m4.resize(3,2);
std::cout << "after resize m5 matrix" << std::endl << m4 << std::endl;    

}