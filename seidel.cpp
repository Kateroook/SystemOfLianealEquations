#include "functions.h"

void seidel(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    vector<double> x(n, 0.0);
    vector<double> x_prev(n, 0.0);
    int iteration = 0;

    printSolution(x,'x',  "Initial approximation");

    double epsilon;
    cout << "Enter precision (epsilon): ";
    cin >> epsilon;

    bool isOk = true;
    while (isOk == true) {
        iteration++;
        x_prev = x;

        for (int i = 0; i < n; i++) {
            double sum1 = 0.0, sum2 = 0.0;

            for (int j = 0; j < i; j++) {
                double prev_sum = sum1;
                sum1 += (A[i][j] /A[i][i])* x[j];
                cout << "sum1^k+1[" << j+1 << "] = sum1^k[" << j+1 << "] + a[" << i+1 << "][" <<j+1 << "] / a[" << i+1 << "][" << i+1 << "] * x^k+1["
                <<j+1 << "] = "<< prev_sum << " + (" << A[i][j] << " / " << A[i][i] << ") * " << x[j] << " = " << sum1 << endl;
            }
            for (int j = i + 1; j < n; j++) {
                double prev_sum = sum2;
                sum2 += (A[i][j]/A[i][i]) * x_prev[j];
                cout << "sum2^k+1[" << i+1 << "] = sum1^k[" << j+1 << "] + a[" << i+1 << "][" <<j+1 << "] / a[" << i+1 << "][" << i+1 << "] * x^k["
                <<j+1 << "] = " << prev_sum << " + (" << A[i][j] << " / " << A[i][i] << ") * " << x_prev[j] << " = " << sum2 << endl;
            }

            x[i] = - sum1 - sum2 + (b[i] / A[i][i]);
            cout << "x[" << i+1 <<"] = - sum1 - sum2 + b[" << i+1 << "]/ a[" << i+1 << "][" << i+1 << "] = "
            << -sum1 << " - " << sum2 << " + " << b[i] << " / " << A[i][i] << " = " << -sum1-sum2 << " + " << b[i] / A[i][i] << " = " << x[i] << endl << endl;
        }

        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            norm += pow(x[i] - x_prev[i], 2);
        }
        norm = sqrt(norm);

        printSolution(x,'x', "Iteration " + to_string(iteration));
        cout << "Stopping condition: ||x^(n) - x^(n-1)|| = " << norm;

        if (norm <= epsilon) {
            cout << " <= " << epsilon << endl << endl;
            isOk = false;
        }
        else {
            cout << " > " << epsilon << endl<< endl;
        }
    }
    printSolution(x, 'x', "Solution after " + to_string(iteration) + " iterations");
}