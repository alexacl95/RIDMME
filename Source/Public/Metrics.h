#ifndef METRIC_H
#define METRIC_H

#include "MatrixRIDMME.h"
#include <cmath>

double MetricEval(Array ydata, Array Constraint, int ErrorMetric);
Array TimeSeriesError(Array (*fnc)(Array X_inf, Array X_sup), Matrix X_inf,
                      Matrix X_sup, Matrix Constraint, int ErrorMetric);
void MultipleError(Array (*fnc)(Array X_inf, Array X_sup), MatrixList X_inf,
                   MatrixList X_sup, Matrix Constraint, int ErrorMetric,
                   string Name);

#endif