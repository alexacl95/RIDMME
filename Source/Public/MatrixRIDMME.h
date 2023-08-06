#ifndef MA_H
#define MA_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

using Array = vector<double>;
using Matrix = std::vector<Array>;
using MatrixList = vector<Matrix>;

using ArrayInt = vector<int>;
using MatrixInt = vector<ArrayInt>;

using ArrayBool = vector<bool>;

void PrintMatrix(Matrix &vec);
void PrintMatrixInt(MatrixInt &vec);
void PrintArray(Array &vec);
void MergeTwoArray(Array X_inf, Array X_sup, Matrix &X);
ArrayInt getCumulativeSum(ArrayInt &arr);
Matrix combineMatrices(const Matrix &Matrix1, const Matrix &Matrix2,
                       const int MinMax);
void AppTwoArray(Array X1, Array x2, Array &combinedArray);
Array sumVectors(Array &Array1, Array &Array2);
Array susVectors(Array &Array1, Array &Array2);
Array MatrixVectorMultiply(const Matrix &Matrix, const Array &Array);
void saveMatrixToCSV(const Matrix &Matrix, const string &filename);
void saveMatricesToCSV(const Matrix &Matrix1, const Matrix &Matrix2,
                       const string &filename);
void saveArrayToCSV(const Array &Array1, const Array &Array2,
                    const string &filename);
bool hasNaN(const Matrix &Matrix);
double *vectorToPrimVector(const Array &vec);
ArrayBool checkForNaN(const Matrix &Matrix1, const Matrix &Matrix2,
                      int RoworCol);
Matrix findMinMaxMatrix(const MatrixList &matrices, int MinMax);
Matrix merge2Vectors(const Array &vector1, const Array &vector2);
#endif