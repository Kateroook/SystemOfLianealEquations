#include "functions.h"

Matrix calculateInverseMatrix(Matrix A) {
    int n = A.size();
    Matrix inverse = createIdentityMatrix(n);

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;
        double maxVal = fabs(A[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(A[i][k]) > maxVal) {
                maxVal = fabs(A[i][k]);
                pivotRow = i;
            }
        }

        if (pivotRow != k) {
            swap(A[k], A[pivotRow]);
            swap(inverse[k], inverse[pivotRow]);
        }

        double pivotValue = A[k][k];
        for (int j = 0; j < n; ++j) {
            A[k][j] /= pivotValue;
            inverse[k][j] /= pivotValue;
            if (std::signbit(inverse[k][j])) {
                if (inverse[k][j] == 0) {
                    inverse[k][j] = 0;
                }
            }
        }
        // printMatrix(A, "A after normalization of pivot row");
        // printMatrix(inverse, "Inverse matrix after elimination of the column");

        for (int i = 0; i < n; ++i) {
            if (i != k) {
                double factor = A[i][k];
                for (int j = 0; j < n; ++j) {
                    A[i][j] -= factor * A[k][j];
                    inverse[i][j] -= factor * inverse[k][j];
                }
            }
        }
        // printMatrix(A, "A after elimination of the column");
        // printMatrix(inverse, "Inverse matrix after elimination of the column");
    }
    return inverse;
}

void gauss(Matrix A, vector<double> b) {
    int n = A.size();
    Matrix P_total = createIdentityMatrix(n);
    Matrix I = createIdentityMatrix(n);
    Matrix A_init = A;
    double det = 1.0;
    int swaps = 0;

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;
        double maxVal = A[k][k];
        for (int i = k + 1; i < n; ++i) {
            if (fabs(A[i][k]) > fabs(maxVal)) {
                maxVal = A[i][k];
                pivotRow = i;
            }
        }
        cout << "Pivot element: " << maxVal << endl << endl;

        // P_k
        Matrix P_k = createIdentityMatrix(n);
        if (pivotRow != k) {
            swap(P_k[k], P_k[pivotRow]);
            swap(b[k], b[pivotRow]);
            ++swaps;
        }
        A = multiplyMatrices(P_k, A);

        printMatrix(P_k, "P[" + to_string(k + 1) + "]");
        P_total = multiplyMatrices(P_k, P_total);

        printMatrix(A,b,"A'[" + to_string(k+1) +"] = P[" + to_string(k+1) + "]* A'[" + to_string(k) + "]");

        //M_k
        Matrix M_k = createIdentityMatrix(n);
        M_k[k][k] = 1 / A[k][k];
        for (int i = k + 1; i < n; ++i) {
            M_k[i][k] = -A[i][k] / A[k][k];
        }
        det *= A[k][k];

        A = multiplyMatrices(M_k, A);
        b = multiplyMatrixVector(M_k, b);

        printMatrix(M_k, "M[" + to_string(k + 1) + "]");
        printMatrix(A,b, "A[" + to_string(k + 1) + "] = M[" + to_string(k+1) + "] * A'[" +to_string(k+1) +"]");
    }
    //printMatrix(P_total, "P");

    // Back substitution
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            cout << "x[" << i+1 << "] = x[" << i+1 << "] - a[" << i+1 << "][" <<j+1 << "] * x[" << j+1 <<"] = " << x[i] << " - " <<A[i][j] << "*" << x[j] << " = " << x[i] << " - " << A[i][j] * x[j];
            x[i] -= A[i][j] * x[j];
            cout << " = " << x[i]<< endl;
        }
        x[i] /= A[i][i];

        cout << "x[" << i+1 << "] = x[" << i+1 << "] / a[" << i+1 << "][" << i+1 << "] = " << x[i] << " / " << A[i][i] << " = " << x[i] << endl;
    }

    I = calculateInverseMatrix(A_init);

    printMatrix(A_init, "A");

    printSolution(x, 'x', "Solution");
    cout << "\nChecking if solution is right: A*x = b" << endl;
    auto solved = multiplyMatrixVector(A_init, x);
    printSolution(solved, 'b', "b");

    cout << endl;
    printMatrix(I, "A^-1");
    cout << "\ndet A = (-1)^p a[i][i] \np = " << swaps << endl;
    cout << "Determinant of A: (-1)^" << swaps << " * "<< det;
    det = (swaps % 2 == 0) ? det : -det;
    cout << " = " << det << endl;
}