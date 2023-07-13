#ifndef MonCar_H
#define MonCar_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "Matrix.h"
#include "AutomaticDifferentiation.h"

double RandN(double low, double up);
void MontecarloSets(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const),
                                   const array& x_low, const array& x_up,
                                   const array& p_low, const array& p_up, int sample_succ_number,
                                   int& time, matrix& Result_inf, matrix& Result_sup);

#endif