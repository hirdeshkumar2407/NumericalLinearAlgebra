#include <iostream>
#include <Eigen/Dense>

 
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
int main(){
MatrixXd m(3,3);

m << 1,2,3,
     4,5,6,
     7,8,9;
/*for (int i=0 ; i<3 ; i++){
    for (int j=0 ; j<3 ; j++){
        if (i==j){
        std:: cout << "Element number " << i+j+1 << "=" << m(i,j) << std::endl;}
        else if (i==1){
            std:: cout << "Element number " << i+j+2 << "=" << m(i,j) << std::endl;
        }
        else if (i==2){
            std:: cout << "Element number " << i+j+3 << "=" << m(i,j) << std::endl;
        }
    }

}*/
std::cout << "m matrix" << std::endl << m << std::endl;
std::cout << "m size ="  << m.size() << std::endl;
std::cout << "m rows ="  << m.rows() << std::endl;
std::cout << "m cols ="  << m.cols() << std::endl;
}