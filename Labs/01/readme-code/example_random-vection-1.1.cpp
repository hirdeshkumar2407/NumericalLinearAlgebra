
#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

int main()
{
    MatrixXd m = MatrixXd::Random(3,3);
    m = (  m + MatrixXd::Constant(3,3,1.0)  ) * 10;

    cout << "m =" << endl << m << endl;
    VectorXd v(3);
    v << 1, 0,0;
    cout << "m * v =" << endl << m * v << endl;
}