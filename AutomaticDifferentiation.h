#ifndef AD_H
#define AD_H

#include "AutoDiff/IntervalAutoDiff.h"
#include "RIDMME_Interval_Autodiff.h"
#include "Matrix.h"
#include <limits>
#include <cmath>
#include <algorithm>

void AutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&,INTERVAL_VECTOR&, int const),
                    array X_inf, array X_sup, array P_inf, array P_sup, matrix& GrdInf, matrix& GrdSup, matrix& GrdInfP, matrix& GrdSupP);
void MultipleAutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const), matrix X_sub, matrixInt allPerms,
                            arrayInt N, matrix& Jx_inf, matrix& Jx_sup, matrix& Jp_inf, matrix& Jp_sup, arrayBool& AnyNanElement);
void ToDoubleEvaluationModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const), double t, array& x, array& p);
RIDMME_INTERVAL_AUTODIFF MergeTwoIntervalArr(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P);
#endif