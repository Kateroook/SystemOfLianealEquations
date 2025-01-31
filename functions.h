#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

using namespace std;

typedef vector<vector<double>> Matrix;

void printMatrix(const Matrix &mat, const string &name);
void printMatrix(const vector<vector<double>>& A, const vector<double>& f, const string &name);
void printVector(const vector<double>& vec, const string &name);
void printSolution(const vector<double>& vec,const char &symbol, const string &name);

Matrix createIdentityMatrix(int size);
Matrix multiplyMatrices(const Matrix &A, const Matrix &B);
Matrix calculateInverseMatrix(Matrix A);
vector<double> multiplyMatrixVector(const vector<vector<double>>& A, const vector<double>& v);
Matrix transposeMatrix(const vector<vector<double>>& A);


void gauss(Matrix A, vector<double> b);
void sqrtMethod(const vector<vector<double>>& A, const vector<double>& b);
void seidel(const vector<vector<double>>& A, const vector<double>& b);

#endif //FUNCTIONS_H