#ifndef AD_H
#define AD_H

#include "AutoDiff/IntervalAutoDiff.h"
#include "MatrixRIDMME.h"
#include "RIDMME_Interval_Autodiff.h"
#include <algorithm>
#include <cmath>
#include <limits>

void AutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                   INTERVAL_VECTOR &,
                                                   int const, int const),
                   Array X_inf, Array X_sup, Array P_inf, Array P_sup,
                   int iterationTime, Matrix &GrdInf, Matrix &GrdSup, Matrix &GrdInfP,
                   Matrix &GrdSupP);
void MultipleAutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                           INTERVAL_VECTOR &,
                                                           int const, int const),
                           Matrix X_sub, int iterationTime, MatrixInt allPerms, ArrayInt N,
                           Matrix &Jx_inf, Matrix &Jx_sup, Matrix &Jp_inf,
                           Matrix &Jp_sup, ArrayBool &AnyNanElement);
void ToDoubleEvaluationModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                             INTERVAL_VECTOR &,
                                                             int const, int const),
                             double t, Array &x, Array &p);
RIDMME_INTERVAL_AUTODIFF MergeTwoIntervalArr(INTERVAL_VECTOR &v,
                                             INTERVAL_VECTOR &P);
#endif