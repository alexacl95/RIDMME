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

void AppTwoArray(array x1, array x2, array& combinedArray){
    combinedArray = x1;
    if(x2.size()!=0){
    combinedArray.insert(combinedArray.end(), x2.begin(), x2.end());
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

matrix combineMatrices(const matrix& matrix1, const matrix& matrix2, const int MinMax){
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

array sumVectors(array& array1, array& array2){
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

array susVectors(array& array1, array& array2){
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

array matrixVectorMultiply(const matrix& Matrix, const array& Array){
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

void saveMatrixToCSV(const matrix& Matrix, const string& filename){
    ofstream file(filename);

    if (!file) {
        cout << "Error: Failed to create the file." << endl;
        return;
    }
    // Iterate over the rows of the matrix
    for (const auto& row : Matrix) {
        // Iterate over the elements of each row
        for (const auto& element : row) {
            // Write the element to the file
            file << element << ",";
        }
        // Move to the next line
        file << endl;
    }

    file.close();
    cout << "Matrix saved to " << filename << endl;
}

void saveMatricesToCSV(const matrix& matrix1, const matrix& matrix2, const string& filename){
    ofstream file(filename);

    // Write matrix1 to CSV
    for (const auto& row : matrix1) {
        for (const auto& element : row) {
            file << element << ",";
        }
        file << endl;
    }
    // Write matrix2 to CSV
    for (const auto& row : matrix2) {
        for (const auto& element : row) {
            file << element << ",";
        }
        file << endl;
    }

    file.close();
}

bool hasNaN(const matrix& Matrix){
    for (const auto& row : Matrix) {
        for (const auto& element : row) {
            if (isnan(element)) {
                return true;
            }
        }
    }
    return false;
}

double* vectorToPrimVector(const array& vec) {
    double* arr = new double[vec.size()];
    std::copy(vec.begin(), vec.end(), arr);
    return arr;
}

arrayBool checkColumnsForNaN(const matrix& Matrix, const matrix& Matrix2){
    arrayBool result(Matrix[0].size(), false);
    for (size_t col = 0; col < Matrix[0].size(); ++col) {
        for (size_t row = 0; row < Matrix.size(); ++row) {
            if (isnan(Matrix[row][col])||isnan(Matrix2[row][col])) {
                result[col] = true;
            }
        }
    }
    return result;
}

matrix findMinMaxMatrix(const matrixList& matrices, int MinMax){
    if (matrices.empty()) {
        return matrix(); // Return an empty matrix if the list is empty
    }    
    const size_t numRows = matrices[0].size();
    const size_t numCols = matrices[0][0].size();
    
    matrix minmaxMatrix(numRows, array(numCols, numeric_limits<double>::infinity()));
    
    for (const auto& matrix : matrices) {
        for (size_t row = 0; row < numRows; ++row) {
            for (size_t col = 0; col < numCols; ++col) {
                if(MinMax == 0){
                    minmaxMatrix[row][col] = min(minmaxMatrix[row][col], matrix[row][col]);
                }else if (MinMax == 1){
                    minmaxMatrix[row][col] = max(minmaxMatrix[row][col], matrix[row][col]);
                }
                
            }
        }
    }    
    return minmaxMatrix;
}