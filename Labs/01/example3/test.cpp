#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
int main()
{
    // Initialize a 100x100 matrix
    MatrixXd m(100, 100);

    // Set all elements of the matrix to zero
    m.setZero();

    // Set the main diagonal (diagonal with offset 0) to 2
   
    
    for(int i=0;i<100;i++)
    {
       m(i,i)=2;

       if (i<99)
       {
           m(i,i+1)=-1;
           m(i+1,i)=-1;
       }
       
   
     
    }



    // Output the matrix
    std::cout << "m = \n" << m << std::endl;
    std::cout << "m.norm  = " << m.norm() << std::endl;
    std::cout << "norm of symmetric part  = " << (m.transpose() +m).norm() << std::endl;
    std::cout << "m bottomRightConner  = " << m.bottomRightCorner(50,50) << std::endl;
   
     VectorXd  v = VectorXd::Constant(50,1.0);
   std::cout << "vector v and 50 A=" << v.dot(m.row(0).head(50)) << std::endl;

    return 0;
}
