#ifndef AD_H
#define AD_H

#include "AutoDiff/IntervalAutoDiff.h"
#include "AutoDiff/AutoDiff.h"
#include "Matrix.h"

INTERVAL_AUTODIFF ModelFunction (INTERVAL_VECTOR & v, int const i);
void AutoDiffModel(array X_inf, array X_sup, matrix& GrdInf, matrix& GrdSup);

#endif