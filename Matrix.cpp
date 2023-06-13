#include "Matrix.h"

void PrintMatrix(matrix& vec){
    for (int i = 0; i < vec.size(); i++){
        for (int j = 0; j < vec[i].size(); j++){
            cout << vec[i][j]<< " ";
        }
        cout << endl;
    }
    cout << endl;
}

void PrintMatrixInt(matrixInt& vec){
    for (int i = 0; i < vec.size(); i++){
        for (int j = 0; j < vec[i].size(); j++){
            cout << vec[i][j]<< " ";
        }
        cout << endl;
    }
    cout << endl;
}

void PrintArray(array& vec){
    for (int i = 0; i < vec.size(); i++){
        cout << vec[i]<< " ";
    }
    cout << endl;
}

void MergeTwoArray(array X_inf, array X_sup, matrix& X){
    int n_x = X_inf.size();
    X.resize(n_x);
    for(int i = 0; i < n_x; i++){
        X[i].resize(2);
        X[i][0] = X_inf[i];
        X[i][1] = X_sup[i];  
    }
}

arrayInt getCumulativeSum(arrayInt &arr){
   arrayInt cumulativeSum{0};
   for(int i = 0; i < arr.size(); i++) {
       int prefixSum = 0;
       for (int j = 0; j <= i; j++) {
           prefixSum += arr[j];
       }
       cumulativeSum.push_back(prefixSum);
   }
   return cumulativeSum;
}

matrix combineMatrices(const matrix& matrix1, const matrix& matrix2, const int MinMax ) {
    double rows = matrix1.size();
    double cols = matrix1[0].size();

    matrix combinedMatrix(rows, array(cols));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if(MinMax == 0){
                combinedMatrix[i][j] = min(matrix1[i][j], matrix2[i][j]);
            }else if (MinMax == 1) {
                combinedMatrix[i][j] = max(matrix1[i][j], matrix2[i][j]);
            }            
        }
    }
    return combinedMatrix;
}

array sumVectors(array& array1, array& array2) {
    // Check if the input vectors have the same size
    if (array1.size() != array2.size()) {
        cout << "Error: Vectors must have the same size." << endl;
        return {};
    }
    array result(array1.size());
    // Iterate over the elements of the vectors and add them together
    for (int i = 0; i < array1.size(); i++) {
        result[i] = array1[i] + array2[i];
    }
    return result;
}

array susVectors(array& array1, array& array2) {
    // Check if the input vectors have the same size
    if (array1.size() != array2.size()) {
        cout << "Error: Vectors must have the same size." << endl;
        return {};
    }
    array result(array1.size());
    // Iterate over the elements of the vectors and add them together
    for (int i = 0; i < array1.size(); i++) {
        result[i] = array1[i] - array2[i];
    }
    return result;
}

array matrixVectorMultiply(const matrix& Matrix, const array& Array) {
    // Check if the matrix and vector dimensions are compatible
    int numRows = Matrix.size();
    int numCols = Matrix[0].size();
    int ArraySize = Array.size();

    if (numCols != ArraySize) {
        cout << "Error: Matrix and vector dimensions are incompatible." << endl;
        return {};
    }
    array result(numRows, 0.0);
    // Perform matrix-vector multiplication
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i] += Matrix[i][j] * Array[j];
        }
    }
    return result;
}
