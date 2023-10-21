#include "../include/matrix.h"

using namespace std;

int main()
{
   
    Matrix<int> A(3, 3, 0);
    Matrix<int> B(2, 2, 0);
    A(0, 0) = 1;
    A(0, 1) = 3;
    A(0, 2) = 9;
    A(1, 0) = 6;
    A(1, 1) = 2;
    A(1, 2) = 3;
    A(2, 0) = 9;
    A(2, 1) = 8;
    A(2, 2) = 5;
    
   
    
    Matrix<double> C;
    C = (A / 2);
    Matrix<double> D = C * 3;

    D.Print();

    const int kN = 3;
    double x[kN] = {2, 3, 4};

    double* b = solveLinearSystem(A, x, kN);
    cout << '(';
    for (int i = 0; i < kN; ++i) {
        cout << '(' << b[i] << ', ';
    }
    cout << ')';

    return 0;
}
