### 2.3 Load and export sparse matrices

The matrix market format is a common way to store spare matrices. It is suported by most of the existing Numerical Linear Algebra libraries (Eigen and LIS in particular). 

To export your matrices and right-hand-side vectors in the matrix-market format, we can use the `unsupported SparseExtra` module.

```
#include <unsupported/Eigen/SparseExtra>
...
Eigen::saveMarket(A, "filename.mtx");
Eigen::saveMarketVector(B, "filename_b.mtx");
```

To load a matrix in the matrix market format, follow the instructions below:

- in the terminal, use `wget` to download a matrix from the matrix market, e.g. `wget https://math.nist.gov/pub/MatrixMarket2/NEP/mhd/mhd416a.mtx.gz`.

- unzip the file by typing `gzip -dk mhd416a.mtx.gz`

- in Eigen, include the `unsupported SparseExtra` module and use 

```
SparseMatrix<double> mat;
loadMarket(mat, "mhd416a.mtx");
```

To export a matrix in the matrix market format, follow the instructions below:

```
std::string matrixFileOut("./matrixName.mtx");
Eigen::saveMarket(mat, matrixFileOut);
```

### Exercise

Compile and test the following example 

```
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    // Load matrix
    SparseMatrix<double> mat;
    loadMarket(mat, "mhd416a.mtx");

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);      
    // define exact solution
    VectorXd b = mat*xe;                 // compute right-hand side
    cout << b << endl;
    
    // Export vector in .mtx format
    int n = b.size();
    // Eigen::saveMarketVector(b, "./rhs.mtx");
    FILE* out = fopen("rhs.mtx","w");
    fprintf(out,"%%%%MatrixMarket vector coordinate real general\n");
    fprintf(out,"%d\n", n);
    for (int i=0; i<n; i++) {
        fprintf(out,"%d %f\n", i ,b(i));
    }
    fclose(out);

    return 0;    
}
```