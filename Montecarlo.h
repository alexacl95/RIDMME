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
                                   int& time, const std::vector<double>& x_low, const std::vector<double>& x_up,
                                   const std::vector<double>& p_low, const std::vector<double>& p_up, int sample_succ_number);

#endif