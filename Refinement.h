#ifndef RE_H
#define RE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Matrix.h"

matrix NonUniformSubdivision(array Factor_inf, array Factor_sup, arrayInt& N);
matrixInt generateAllPermutations(arrayInt& N);
void ReduceN(arrayInt& N, arrayInt& N_0, int TimesN0, double PercentN);
void refinementIteration();

#endif