#include "MatrixRIDMME.h"

void PrintMatrix(Matrix &vec) {
  for (int i = 0; i < vec.size(); i++) {
    for (int j = 0; j < vec[i].size(); j++) {
      cout << vec[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void PrintMatrixInt(MatrixInt &vec) {
  for (int i = 0; i < vec.size(); i++) {
    for (int j = 0; j < vec[i].size(); j++) {
      cout << vec[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void PrintArray(Array &vec) {
  for (int i = 0; i < vec.size(); i++) {
    cout << vec[i] << " ";
  }
  cout << endl;
}

void MergeTwoArray(Array X_inf, Array X_sup, Matrix &X) {
  int n_x = X_inf.size();
  X.resize(n_x);
  for (int i = 0; i < n_x; i++) {
    X[i].resize(2);
    X[i][0] = X_inf[i];
    X[i][1] = X_sup[i];
  }
}

void AppTwoArray(Array x1, Array x2, Array &combinedArray) {
  combinedArray = x1;
  if (x2.size() != 0) {
    combinedArray.insert(combinedArray.end(), x2.begin(), x2.end());
  }
}

ArrayInt getCumulativeSum(ArrayInt &arr) {
  ArrayInt cumulativeSum{0};
  for (int i = 0; i < arr.size(); i++) {
    int prefixSum = 0;
    for (int j = 0; j <= i; j++) {
      prefixSum += arr[j];
    }
    cumulativeSum.push_back(prefixSum);
  }
  return cumulativeSum;
}

Matrix combineMatrices(const Matrix &Matrix1, const Matrix &Matrix2,
                       const int MinMax) {
  double rows = Matrix1.size();
  double cols = Matrix1[0].size();

  Matrix combinedMatrix(rows, Array(cols));

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (MinMax == 0) {
        combinedMatrix[i][j] = min(Matrix1[i][j], Matrix2[i][j]);
      } else if (MinMax == 1) {
        combinedMatrix[i][j] = max(Matrix1[i][j], Matrix2[i][j]);
      }
    }
  }
  return combinedMatrix;
}

Array sumVectors(Array &Array1, Array &Array2) {
  // Check if the input vectors have the same size
  if (Array1.size() != Array2.size()) {
    cout << "Error: Vectors must have the same size." << endl;
    return {};
  }
  Array result(Array1.size());
  // Iterate over the elements of the vectors and add them together
  for (int i = 0; i < Array1.size(); i++) {
    result[i] = Array1[i] + Array2[i];
  }
  return result;
}

Array susVectors(Array &Array1, Array &Array2) {
  // Check if the input vectors have the same size
  if (Array1.size() != Array2.size()) {
    cout << "Error: Vectors must have the same size." << endl;
    return {};
  }
  Array result(Array1.size());
  // Iterate over the elements of the vectors and add them together
  for (int i = 0; i < Array1.size(); i++) {
    result[i] = Array1[i] - Array2[i];
  }
  return result;
}

Array MatrixVectorMultiply(const Matrix &MatrixInput, const Array &ArrayInput) {
  // Check if the Matrix and vector dimensions are compatible
  int numRows = MatrixInput.size();
  int numCols = MatrixInput[0].size();
  int ArraySize = ArrayInput.size();

  if (numCols != ArraySize) {
    cout << "Error: Matrix and vector dimensions are incompatible." << endl;
    return {};
  }
  Array result(numRows, 0);
  // Perform Matrix-vector multiplication
  for (int i = 0; i < numRows; i++) {
    for (int j = 0; j < numCols; j++) {
      result[i] += MatrixInput[i][j] * ArrayInput[j];
    }
  }
  return result;
}

void saveMatrixToCSV(const Matrix &Matrix, const string &filename) {
  ofstream file(filename);

  if (!file) {
    cout << "Error: Failed to create the file." << endl;
    return;
  }
  // Iterate over the rows of the Matrix
  for (const auto &row : Matrix) {
    // Iterate over the elements of each row
    for (const auto &element : row) {
      // Write the element to the file
      file << element << ",";
    }
    // Move to the next line
    file << endl;
  }

  file.close();
  cout << "Matrix saved to " << filename << endl;
}

void saveMatricesToCSV(const Matrix &Matrix1, const Matrix &Matrix2,
                       const string &filename) {
  ofstream file(filename);

  // Write Matrix1 to CSV
  for (const auto &row : Matrix1) {
    for (const auto &element : row) {
      file << element << ",";
    }
    file << endl;
  }
  // Write Matrix2 to CSV
  for (const auto &row : Matrix2) {
    for (const auto &element : row) {
      file << element << ",";
    }
    file << endl;
  }

  file.close();
}
void saveArrayToCSV(const Array &Array1, const Array &Array2,
                    const string &filename) {
  ofstream file(filename);

  // Write Matrix1 to CSV
  for (const auto &row : Array1) {
    file << row << "," << endl;
  }
  // Write Matrix2 to CSV
  for (const auto &row : Array2) {
    file << row << "," << endl;
  }

  file.close();
}

bool hasNaN(const Matrix &Matrix) {
  for (const auto &row : Matrix) {
    for (const auto &element : row) {
      if (isnan(element)) {
        return true;
      }
    }
  }
  return false;
}

double *vectorToPrimVector(const Array &vec) {
  double *arr = new double[vec.size()];
  std::copy(vec.begin(), vec.end(), arr);
  return arr;
}

ArrayBool checkForNaN(const Matrix &Matrix1, const Matrix &Matrix2,
                      int RoworCol) {
  ArrayBool result(Matrix1[0].size(), false);
  for (size_t col = 0; col < Matrix1[0].size(); ++col) {
    for (size_t row = 0; row < Matrix1.size(); ++row) {
      if (isnan(Matrix1[row][col]) || isnan(Matrix2[row][col])) {
        if (RoworCol == 1) {
          result[row] = true;
        } else if (RoworCol == 2) {
          result[col] = true;
        }
      }
    }
  }
  return result;
}

Matrix findMinMaxMatrix(const MatrixList &matrices, int MinMax) {
  if (matrices.empty()) {
    return Matrix(); // Return an empty Matrix if the list is empty
  }
  const size_t numRows = matrices[0].size();
  const size_t numCols = matrices[0][0].size();

  Matrix minmaxMatrix;
  if (MinMax == 0) {
    minmaxMatrix =
        Matrix(numRows, Array(numCols, numeric_limits<double>::infinity()));
  } else if (MinMax == 1) {
    minmaxMatrix =
        Matrix(numRows, Array(numCols, -numeric_limits<double>::infinity()));
  }

  for (const auto &Matrix : matrices) {
    for (size_t row = 0; row < numRows; ++row) {
      for (size_t col = 0; col < numCols; ++col) {
        if (MinMax == 0) {
          minmaxMatrix[row][col] =
              min(minmaxMatrix[row][col], Matrix[row][col]);
        } else if (MinMax == 1) {
          minmaxMatrix[row][col] =
              max(minmaxMatrix[row][col], Matrix[row][col]);
        }
      }
    }
  }
  return minmaxMatrix;
}

Matrix merge2Vectors(const Array &vector1, const Array &vector2) {
  size_t size = vector1.size();
  if (size != vector2.size()) {
    cout << "Error: Vectors have different sizes." << endl;
    return Matrix();
  }

  Matrix result(size, Array(2));
  for (size_t i = 0; i < size; ++i) {
    result[i][0] = vector1[i];
    result[i][1] = vector2[i];
  }

  return result;
}