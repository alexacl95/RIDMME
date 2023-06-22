#ifndef MA_H
#define MA_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

using array = vector<double>;
using matrix = vector<array>;

using arrayInt = vector<int>;
using matrixInt = vector<arrayInt>;

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
bool hasNaN(const matrix& Matrix);
double* vectorToPrimVector(const array& vec);

#endif