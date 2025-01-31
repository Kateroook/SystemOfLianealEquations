#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

bool isSymmetric(const std::vector<std::vector<double>>& matrix) {
    size_t n = matrix.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}
bool hasDiagonalDominance(const std::vector<std::vector<double>>& matrix) {
    size_t n = matrix.size();
    for (size_t i = 0; i < n; ++i) {
        double diagonal = std::abs(matrix[i][i]);
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += std::abs(matrix[i][j]);
            }
        }
        if (diagonal < sum) {
            return false;
        }
    }
    return true;
}
void checkMatrixProperties(const std::vector<std::vector<double>>& matrix) {
    std::cout << "Checking matrix properties...\n";

    if (isSymmetric(matrix)) {
        std::cout << "The matrix is symmetric.\n";
    } else {
        std::cout << "The matrix is not symmetric.\n";
    }

    if (hasDiagonalDominance(matrix)) {
        std::cout << "The matrix has diagonal dominance.\n";
    } else {
        std::cout << "The matrix does not have diagonal dominance.\n";
    }
}

void printMatrix(const vector<vector<double>>& A, const vector<double>& f) {
int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << A[i][j] << " ";
        }
        cout << "| " << setw(10) << f[i] << endl;
    }
    cout << endl;
}

void printMatrix(const vector<vector<double>>& A) {
    for (const auto& row : A) {
        for (double val : row) {
            cout << setw(10) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void printVector(const vector<double>& vec) {
    for (double val : vec) {
        cout << setw(10) << val << " ";
    }
    cout << endl;
}

double calculateDeterminant(vector<vector<double>> A, int n) {
    double det = 1.0;

    for (int i = 0; i < n; i++) {
        double maxEl = abs(A[i][i]);
        int row = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                row = k;
            }
        }

        if (row != i) {
            swap(A[i], A[row]);
            det = -det;
        }

        det *= A[i][i];

        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }

    return det;
}

vector<vector<double>> inverseMatrix(vector<vector<double>> A, int n) {
    vector<vector<double>> inv(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        inv[i][i] = 1.0;
    }

    for (int i = 0; i < n; i++) {
        double maxEl = abs(A[i][i]);
        int row = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                row = k;
            }
        }
        if (row != i) {
            swap(A[i], A[row]);
            swap(inv[i], inv[row]);
        }

        double pivot = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            for (int k = 0; k < n; k++) {
                inv[j][k] -= factor * inv[i][k];
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double factor = A[j][i];
            for (int k = 0; k < n; k++) {
                inv[j][k] -= factor * inv[i][k];
            }
        }
    }
    return inv;
}

vector<vector<double>> transposeMatrix(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> T(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            T[j][i] = A[i][j];
        }
    }
    return T;
}

vector<vector<double>> multiplyMatrices(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int n = A.size();
    vector<vector<double>> result(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

vector<double> multiplyMatrixVector(const vector<vector<double>>& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// Function to check the sufficient conditions for convergence
bool checkConvergenceConditions(const vector<vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += abs(A[i][j]);
            }
        }
        if (abs(A[i][i]) < sum) {
            cout << "Row " << i + 1 << ": |a[" << i + 1 << "][" << i + 1 << "]| < Sum of |a[" << i + 1 << "][j]| for j != " << i + 1 << endl;
            return false;
        }
    }
    cout << "Matrix satisfies the sufficient conditions for convergence.\n";
    return true;
}

// Gaussian elimination method to solve the system of equations
vector<double> gaussMethod(vector<vector<double>> A, vector<double> f, int n) {
    vector<double> x(n, 0.0);

    for (int i = 0; i < n; i++) {
        double maxEl = abs(A[i][i]);
        int row = i;
        for (int k = i + 1; k < n; k++) {
            cout << "A[k][i] = A[" << k << "][" << i << "] = " << A[k][i] << endl;
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                row = k;
            }
        }
        std::cout << "Max element in " << i << " row is " << maxEl << std::endl;

        if (row != i) {
            swap(A[i], A[row]);
            swap(f[i], f[row]);
        }

        cout << "Step " << i + 1 << ": After row swap (if needed):" << endl;
        printMatrix(A, f);

        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            f[j] -= factor * f[i];

            cout << "Step " << i + 1 << "." << j + 1 << ": After elimination of row " << j + 1 << ":" << endl;
            printMatrix(A, f);
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = f[i];

        for (int j = i + 1; j < n; j++) {
            cout << "x[" << i+1 << "] = x[" << i+1 << "] - A[" << i+1 << "][" << j+1 << "] * x[" << j+1 << "] = " << x[i] << " - " << A[i][j] << " * " <<  x[j] << " = " << x[i]-A[i][j]*x[j] << endl;
            x[i] -= A[i][j] * x[j];

        }
        cout << "x[" << i+1 << "] = x[" << i+1 << "] / A[" << i+1 << "][" << i+1 << "] = " << x[i] <<" / " << A[i][i] << " = " << x[i] / A[i][i] <<endl;
        x[i] /= A[i][i];
        cout << "Step " << n + n - i << ": After back substitution for x" << i + 1 << ": " << x[i] << endl;
    }
    return x;
}

void seidelMethod(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    vector<double> x(n, 0.0);  // Initial guess (x^0)
    vector<double> x_prev(n, 0.0);
    int iteration = 0;

    double epsilon;
    cout << "Enter precision (epsilon): ";
    cin >> epsilon;

    bool isOk = true;
    while (isOk == true) {
        iteration++;
        x_prev = x;

        // Perform iteration
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0, sum2 = 0.0;

            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x[j];
            }
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x_prev[j];
            }

            x[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        // Check for convergence
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            norm += pow(x[i] - x_prev[i], 2);
        }
        norm = sqrt(norm);

        cout << "Iteration " << iteration << ": ";
        printVector(x);
        cout << "  Stopping condition: ||x^(n) - x^(n-1)|| = " << norm << endl << endl;

        if (norm <= epsilon) {
            isOk = false;
        }
    }

    // Output the solution
    cout << "\nSolution after " << iteration << " iterations:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << "; ";
    }
}

// Square root method
void sqrtMethod(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    vector<vector<double>> S(n, vector<double>(n, 0));
    vector<double> D(n, 0);

    // Decompose A into S and D
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int p = 0; p < i; p++) {
            sum += S[p][i] * S[p][i] * D[p];
        }
        D[i] = (A[i][i] - sum > 0) ? 1 : -1; // Diagonal element of D
        S[i][i] = sqrt(abs(A[i][i] - sum)); // Diagonal element of S

        for (int j = i + 1; j < n; j++) {
            sum = 0.0;
            for (int p = 0; p < i; p++) {
                sum += S[p][i] * D[p] * S[p][j];
            }
            S[i][j] = (A[i][j] - sum) / (D[i] * S[i][i]);
        }
    }

    // Calculate determinant
    double detA = 1.0;
    for (int i = 0; i < n; i++) {
        detA *= D[i] * S[i][i] * S[i][i];
    }
    cout << "Determinant of A: " << detA << endl;

    // Display S and D
    cout << "Matrix S:" << endl;
    printMatrix(S);
    cout << "Matrix D (diagonal elements):" << endl;
    printVector(D);



    // S^T
    vector<vector<double>> ST = transposeMatrix(S);
    cout << "Matrix S^T:" << endl;
    printMatrix(ST);

    // S^T * D
    vector<vector<double>> STD(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            STD[i][j] = ST[i][j] * D[j];
        }
    }
    cout << "Matrix S^T * D:" << endl;
    printMatrix(STD);

    cout << "Matrix S^T * D = f:" << endl;
    printMatrix(STD,b);

    // S^T * D * y = b
    vector<double> y(n, 0);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < i; k++) {
            sum += STD[i][k] * y[k];
        }
        y[i] = (b[i] - sum) / STD[i][i];
    }
    cout << "Vector y:" << endl;
    printVector(y);

    cout << "Matrix S*x = y:" << endl;
    printMatrix(transposeMatrix(STD),y);

    // Solve S * x = y
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int k = i + 1; k < n; k++) {
            sum += S[i][k] * x[k];
        }
        x[i] = (y[i] - sum) / S[i][i];
    }

    // Display solution
    cout << "Solution to the system:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << "; ";
    }
}


int main() {
    // Matrix A (4x4)
    vector<vector<double>> A = {
         {7, -2, 2, -2},
         {-2, 5, -1, 1},
         {2, -1, 6, -2},
         {-2, 1, -2, 6}
    };

    checkMatrixProperties(A);
    // Right-hand side vector b
    vector<double> f = {1, 2, 3, 4};

    // Number of equations (size of the matrix)
    int n = 4;

    cout << "Initial Matrix A:\n";
    printMatrix(A, f);

    // Solve the system using Gauss method
//    vector<double> x = gaussMethod(A, f, n);
//
//    // Calculate the determinant
//    double det = calculateDeterminant(A, n);
//    cout << "\nDeterminant of the matrix: " << det << endl;
//
//    if (det == 0) {
//        cout << "\nThe matrix is singular and does not have an inverse." << endl;
//    }
//    else {
//        // Calculate the inverse of the matrix
//        vector<vector<double>> inv = inverseMatrix(A, n);
//
//        cout << "\nInverse matrix:" << endl;
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                cout << setw(12) << inv[i][j] << " ";
//            }
//            cout << endl;
//        }
//
//
//        // Print the solution
//        cout << "\nSolution to the system:" << endl;
//        for (int i = 0; i < n; i++) {
//            cout << "x" << i + 1 << " = " << x[i] << "; ";
//        }
//    }

    sqrtMethod(A,f);
    //seidelMethod(A,f);
    return 0;
}
