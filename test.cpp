#include <iostream>

#include "functions.h"

int main() {
    Matrix A = {
        {7, -2, 2, -2},
        {-2, 5, -1, 1},
        {2, -1, 6, -2},
        {-2, 1, -2, -2}
    };

    vector<double> b = {1, 2, 3,4};
    printMatrix(A, b, "Initial A");
    //gauss(A, b);
    //sqrtMethod(A, b);
    seidel(A,b);

    return 0;
}