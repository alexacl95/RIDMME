#ifndef METRIC_H
#define METRIC_H

#include <cmath>
#include "Matrix.h"

double MetricEval(array ydata, array Constraint, int ErrorMetric);
array TimeSeriesError(array (*fnc)(array X_inf, array X_sup), array X_inf, array X_sup, array Constraint, int ErrorMetric);

#endif