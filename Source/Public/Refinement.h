#ifndef RE_H
#define RE_H

#include "MatrixRIDMME.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

Matrix NonUniformSubdivision(Array Factor_inf, Array Factor_sup, ArrayInt &N);
MatrixInt generateAllPermutations(ArrayInt &N);
void ReduceN(ArrayInt &N, ArrayInt &N_0, int TimesN0, double PercentN);
void RefinementIteration(ArrayInt &N, int &i, ArrayInt IterationRep, int &Step,
                         int Step0, ArrayBool AnyNanElement, int n_x, int n_p);
#endif