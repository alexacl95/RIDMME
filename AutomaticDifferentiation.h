#ifndef AD_H
#define AD_H

#include "AutoDiff/IntervalAutoDiff.h"
#include "AutoDiff/AutoDiff.h"
#include "Matrix.h"
#include <limits>

void AutoDiffModel(INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), array X_inf, array X_sup, matrix& GrdInf, matrix& GrdSup);
void MultipleAutoDiffModel(INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), matrix X_sub, matrixInt allPerms, arrayInt N, matrix& Jx_inf, matrix& Jx_sup);
void ToDoubleEvaluationModel(INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), double t, array& x, array p);
#endif