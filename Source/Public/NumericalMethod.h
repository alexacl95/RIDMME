#ifndef NM_H
#define NM_H

#include "AutomaticDifferentiation.h"
#include "MatrixRIDMME.h"
#include "Refinement.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

void MixMonotonicityMethod(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                           INTERVAL_VECTOR &,
                                                           int const, int const),
                           double t, Array &x_low, Array &x_up, Array p_low,
                           Array p_up, Matrix J_x_low, Matrix J_x_up,
                           Matrix J_p_low, Matrix J_p_up);
void Simulation(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                INTERVAL_VECTOR &,
                                                int const, int const),
                Array X_inf, Array X_sup, Array P_inf, Array P_sup, ArrayInt N,
                int t_f, Matrix &Results_low, Matrix &Results_up);
int CheckSimSettings(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                     INTERVAL_VECTOR &,
                                                     int const, int const),
                     Array X_inf, Array X_sup, Array P_inf, Array P_sup,
                     ArrayInt N);
void SingleSimulation(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                      INTERVAL_VECTOR &,
                                                      int const, int const),
                      Array &X_inf, Array &X_sup, Array P_inf, Array P_sup,
                      ArrayInt N, int t_f, int i, ArrayBool &AnyNanElement);
void SimulationMixMon(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                      INTERVAL_VECTOR &,
                                                      int const, int const),
                      Array X_inf, Array X_sup, Array P_inf, Array P_sup,
                      int t_f, Matrix &Results_low, Matrix &Results_up);
#endif