#ifndef AD_H
#define AD_H

#include "AutoDiff/IntervalAutoDiff.h"
#include "AutoDiff/AutoDiff.h"
#include "Matrix.h"
#include <limits>

INTERVAL_AUTODIFF ModelFunction (INTERVAL_VECTOR & v, int const i);
void AutoDiffModel(array X_inf, array X_sup, matrix& GrdInf, matrix& GrdSup);
void MultipleAutoDiffModel(matrix X_sub, matrixInt allPerms, arrayInt N, matrix& GrdInf, matrix& GrdSup);
void ToDoubleEvaluationModel(double t, array& x, array p);

#endif