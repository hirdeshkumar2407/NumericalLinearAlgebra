#include <Eigen/Dense>
#include <iostream>
 
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

int main()
{
VectorXd v(6);
v << 1, 2, 3, 4, 5, 6;
cout << "v.head(3) =" << endl << v.head(3) << endl;
cout << "v.tail(3) =" << endl << v.tail(3) << endl;
v.segment(1,4) *= 2;
cout << "v after doubling segment(1,4) =" << endl << v << endl;


/*Eigen provides a set of block operations designed specifically for the special case of vectors:
- Block containing the first $n$ elements: `vector.head(n)`
- Block containing the first $n$ elements: `vector.tail(n)`
- Block containing $n$ elements, starting at position $i$: `vector.segment(i,n)`
*/
/*Eigen also provides special methods for blocks that are flushed against one of the corners or sides of a matrix. For instance:
- Top-left p by q block: `matrix.topLeftCorner(p,q)`
- Bottom-left p by q block: `matrix.bottomLeftCorner(p,q)`
- Top-right p by q block: `matrix.topRightCorner(p,q)`
- Bottom-right p by q block: `matrix.bottomRightCorner(p,q)`

Individual columns and rows are special cases of blocks. Eigen provides methods to easily address them (`.col()` and `.row()`). The argument is the index of the column or row to be accessed. As always in Eigen, indices start at 0.
*/


MatrixXd A = MatrixXd::Random(9,9);
MatrixXd B = A.topLeftCorner(3,6);

VectorXd w = B * v;
cout << "norm of B * v: " << w.norm() << endl;
}


