#ifndef NM_H
#define NM_H

#include "Matrix.h"
#include "AutomaticDifferentiation.h"
#include "Refinement.h"
#include <cmath>
#include <algorithm>
#include <numeric>

void MixMonotonicityMethod(
    INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const),
    double t,
    array& x_low, array& x_up,
    array p_low, array p_up,
    matrix J_x_low, matrix J_x_up,
    matrix J_p_low, matrix J_p_up
);
void Simulation(INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), array X_inf, array X_sup, arrayInt N, int t_f);
#endif