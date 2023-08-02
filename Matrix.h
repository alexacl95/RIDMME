#ifndef MA_H
#define MA_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

using array = vector<double>;
using matrix = vector<array>;
using matrixList = vector<matrix>;

using arrayInt = vector<int>;
using matrixInt = vector<arrayInt>;

using arrayBool = vector<bool>;

void PrintMatrix(matrix& vec);
void PrintMatrixInt(matrixInt& vec);
void PrintArray(array& vec);
void MergeTwoArray(array X_inf, array X_sup, matrix& X);
arrayInt getCumulativeSum(arrayInt &arr);
matrix combineMatrices(const matrix& matrix1, const matrix& matrix2, const int MinMax);
void AppTwoArray(array X1, array x2, array& combinedArray);
array sumVectors(array& array1, array& array2);
array susVectors(array& array1, array& array2);
array matrixVectorMultiply(const matrix& Matrix, const array& Array);
void saveMatrixToCSV(const matrix& Matrix, const string& filename);
void saveMatricesToCSV(const matrix& matrix1, const matrix& matrix2, const string& filename);
void saveArrayToCSV(const array& array1, const array& array2, const string& filename);
bool hasNaN(const matrix& Matrix);
double* vectorToPrimVector(const array& vec);
arrayBool checkForNaN(const matrix& Matrix, const matrix& Matrix2, int RoworCol);
matrix findMinMaxMatrix(const matrixList& matrices, int MinMax);
matrix merge2Vectors(const array& vector1, const array& vector2);
#endif