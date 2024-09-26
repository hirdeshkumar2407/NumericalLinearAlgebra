    #include <Eigen/Dense>
    #include <iostream>
    
    using namespace std;
    using Eigen::VectorXd;
    using Eigen::MatrixXd;

    int main()
    {
    VectorXd v(6);
    cout << "v.head(3) =" << v.head(3) << endl;
    v.segment(1,4) *= 2;
    cout << "v.segment(1,4) =" << v.segment(1,4) << endl;
    MatrixXd A = MatrixXd::Random(9,9);
    MatrixXd B = A.topLeftCorner(3,6);
    VectorXd w = B*v;
    cout << "Matrix A = " << A << endl;
        cout << "Matrix B = " << B << endl;
    cout << "norm of B*v = " << w.norm() << endl;
    }