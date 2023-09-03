#ifndef MonCar_H
#define MonCar_H

#include "AutomaticDifferentiation.h"
#include "MatrixRIDMME.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

double RandN(double low, double up);
void MontecarloSets(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                    INTERVAL_VECTOR &,
                                                    int const, int const),
                    const Array &x_low, const Array &x_up, const Array &p_low,
                    const Array &p_up, int sampleN, int &time,
                    Matrix &Result_inf, Matrix &Result_sup);

#endif