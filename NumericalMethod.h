#ifndef NM_H
#define NM_H

#include "Matrix.h"
#include <cmath>
#include <algorithm>
#include <numeric>

array MonotonicityModel(
    void (*fnc)(double, array&, array),
    double t,
    array x_low, array x_up,
    array p_low, array p_up,
    matrix J_x_low, matrix J_x_up,
    matrix J_p_low, matrix J_p_up
);
#endif