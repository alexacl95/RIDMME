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
void RefinementIteration(array& X_inf, array& X_sup, matrix Results_low, matrix Results_up, arrayInt& N, int& i, arrayInt IterationRep, int& Step, 
                         int Step0, arrayBool AnyNanElement, int n_x, int n_p);                        
#endif