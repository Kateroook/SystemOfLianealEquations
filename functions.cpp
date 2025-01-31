#include "functions.h"

void printMatrix(const Matrix &mat, const string &name) {
    cout << name << " =\n";
    for (const auto &row : mat) {
        for (double val : row) {
            cout << setw(10) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printMatrix(const Matrix& A, const vector<double>& f, const string &name) {
    cout << name << " =\n";
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << A[i][j] << " ";
        }
        cout << setw(5) << "| " << setw(5) << f[i] << endl;
    }
    cout << endl;
}

void printVector(const vector<double>& vec, const string &name) {
    cout << name << "= \n";
    for (double val : vec) {
        cout << setw(10) << val << " ";
    }
    cout << endl;
}

Matrix createIdentityMatrix(int size) {
    Matrix I(size, vector<double>(size, 0.0));
    for (int i = 0; i < size; ++i) {
        I[i][i] = 1.0;
    }
    return I;
}

Matrix multiplyMatrices(const Matrix &A, const Matrix &B) {
    int n = A.size();
    Matrix result(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

vector<double> multiplyMatrixVector(const Matrix& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}


Matrix transposeMatrix(const Matrix& A) {
    int n = A.size();
    Matrix T(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            T[j][i] = A[i][j];
        }
    }
    return T;
}

void printSolution(const vector<double> &vec, const char& symbol, const string &name) {
    cout << name << "= \n";
    for (int i = 0; i < vec.size(); i++) {
        cout <<  symbol << i + 1 << " = " << vec[i] << "; ";
    }
    cout << endl;
}