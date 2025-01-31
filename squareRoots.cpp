#include "functions.h"

void sqrtMethod(const Matrix& A, const vector<double>& b) {
    const int n = A.size();

    cout  << "Checking if matrix is symetric..." << endl;
    const Matrix AT =transposeMatrix(A);
    printMatrix(AT, "A^T");

    for(size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (A[i][j] != AT[i][j]) {
                cout << "Matrix isn't symmetric, can't apply square root method\n";
                return;
            }
        }
    }
    cout << "A = A^T, so we can apply square root method" << endl;


    Matrix S(n, vector<double>(n, 0));
    Matrix D(n, vector<double>(n, 0));


    cout << "\nDecomposition A into S and D: \n";
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int p = 0; p < i; ++p) {
            sum += S[p][i] * S[p][i] * D[p][p];
        }

        D[i][i] = (A[i][i] - sum >= 0) ? 1 : -1;
        cout << "d[" << i+1 << "][" << i+1 <<"] = sgn(a[" << i+1 << "][" << i+1 <<"] - s^2[" << i << "][" << i+1 <<"] * d[" << i << "][" << i
            <<"]) = sgn(" << A[i][i] << " - " <<  sum  << ") = sgn(" << A[i][i] - sum << ") = " << D[i][i] << endl;;

        S[i][i] = sqrt(abs(A[i][i] - sum));
        cout <<  "s[" << i+1 << "][" << i+1 <<"] = sqrt(|a[" << i+1 << "][" << i+1 <<"] - sum s^2[p][i]*d[p][p]|) = sqrt(|" << A[i][i]
            << " - " <<  sum << "|) = sqrt(|" << A[i][i] - sum << "|) = " << S[i][i] << endl;

        for (int j = i + 1; j < n; ++j) {
            sum = 0.0;
            for (int p = 0; p < i; ++p) {
                sum += S[p][i] * D[p][p] * S[p][j];
            }
            cout << "s[" << i+1 << "][" << j+1 <<"] = (a[" << i+1 << "][" << i+1 <<"] - sum(s[p][i]*d[p][p]*s[p][j]) / (D[" << i+1 << "][" << i+1
                <<"] * S[" << i+1 << "][" << i+1 <<"]) = (" <<  A[i][j] << " - " << sum << ") / (" << D[i][i] << " * " << S[i][i] << ") = ";
            S[i][j] = (A[i][j] - sum) / (D[i][i] * S[i][i]);
            cout << S[i][j] << endl;
        }
    }

    // Display S and D
    cout << endl;
    printMatrix(S, "S");
    printMatrix(D, "D");

    // S^T
    vector<vector<double>> ST = transposeMatrix(S);
    printMatrix(ST, "S^T");

    Matrix STD(n, vector<double>(n, 0));
    STD = multiplyMatrices(ST,D);
    printMatrix(STD, "S^T * D");

    printMatrix(STD,b, "S^T * D = b");

    // S^T * D * y = b
    vector<double> y(n, 0);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < i; k++) {
            sum += STD[i][k] * y[k];
        }
        y[i] = (b[i] - sum) / STD[i][i];
        cout << "y[" << i+1 << "] = (b[" << i+1 << "] - sum(STD[i][k] * y[k])) / STD[" << i+1 << "][" << i+1 <<"] = ("
        << b[i] << " - " << sum << ") / " << STD[i][i] << " = " << y[i] << endl;
    }
    //printVector(y,"y");
    printSolution(y, 'y', "y");

    cout << endl;
    printMatrix(transposeMatrix(STD),y, "S*x = y");
    // S * x = y
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int k = i + 1; k < n; k++) {
            sum += S[i][k] * x[k];
        }
        x[i] = (y[i] - sum) / S[i][i];
        cout << "x[" << i+1 << "] = (y[" << i+1 << "] - sum(S[i][k] * x[k])) / S[" << i+1 << "][" << i+1 <<"] = ("
       << y[i] << " - " << sum << ") / " << S[i][i] << " = " << x[i] << endl;
    }
    cout << endl;
    printSolution(x, 'x', "Solution");

    //Calculate determinant
       double detA = 1.0;
        cout << "Determinant of A:\n detA = mult from 1 to " << n << " d[i][i]*s^2[i][i] = ";
       for (int i = 0; i < n; i++) {
           cout << D[i][i] << " * (" << S[i][i] <<")^2 *";
           detA *= D[i][i] * S[i][i] * S[i][i];
       }
       cout << " = " << detA << endl;
}