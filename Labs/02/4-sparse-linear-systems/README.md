### 2.2 Sparse linear systems

<h5>We are focuse on row major, but by default eigen follows column major</h5> 

The class `SparseMatrix` is the main sparse matrix representation of Eigen's sparse module; it offers high performance and low memory usage. It implements a more versatile variant of the widely-used Compressed Column (or Row) Storage scheme.

The `SparseMatrix` and `SparseVector` classes take three template arguments: the scalar type (e.g., double) the storage order (ColMajor or RowMajor, the default is ColMajor) the inner index type (default is int). As for dense Matrix objects, constructors takes the size of the object. Here are some examples:

```
SparseMatrix<double,RowMajor> mat(1000,2000);              
// declares a 1000x2000 row-major compressed sparse matrix of double

SparseVector<double,RowMajor> vec(1000);                   
// declares a row sparse vector of double of size 1000
```

In Eigen, there are several methods available to solve linear systems whenever the matrix is sparse. Because of the special representation of this class of matrices, special care should be taken in order to get a good performance. [This page](https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html) lists the sparse solvers available in Eigen. All the available solvers follow the same general concept.

```
#include <Eigen/RequiredModuleName>
// ...
SparseMatrix<double> A;
// fill A
VectorXd b, x;
// fill b
// solve Ax = b
SolverClassName<SparseMatrix<double> > solver;
solver.compute(A);
if(solver.info()!=Success) {
  // decomposition failed
  return;
}
x = solver.solve(b);
if(solver.info()!=Success) {
  // solving failed
  return;
}
// solve for another right hand side:
x1 = solver.solve(b1);
```
Choseky factorization 
A simple example:



N.B. Block operations works also for sparse matrix, but it is recommended not to use for WRITING or MODIFYING existing matrices. We will only use these operations to extract blocks.
