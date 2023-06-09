#ifndef MA_H
#define MA_H

#include <iostream>
#include <vector>

using namespace std;

using array = vector<double>;
using matrix = vector<array>;

using arrayInt = vector<int>;
using matrixInt = vector<arrayInt>;

void PrintMatrix(matrix& vec);
void PrintMatrixInt(matrixInt& vec);

#endif