
#include "AutomaticDifferentiation.h"
#include "Functions.h"
#include "Interval.h"
#include "IntervalVector.h"
#include "MatrixRIDMME.h"
#include "Metrics.h"
#include "Models.h"
#include "Montecarlo.h"
#include "NumericalMethod.h"
#include "Refinement.h"
#include <iostream>
#include <stdio.h>
#include <vector>

// Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)
// Note: config/systemFolderConfi/BiasRnd.h comment line 42 to compute Nan-inf.

Array SumModel(Array X_inf, Array X_sup) {
  int Size = X_inf.size();
  Array SumStates(2, 0);
  for (int i = 0; i < Size; i++) {
    SumStates[0] = X_inf[i] + SumStates[0];
    SumStates[1] = X_sup[i] + SumStates[1];
  }
  return SumStates;
}

Array SumModelDengue(Array X_inf, Array X_sup) {
  int Size = 3;
  Array SumStates(2, 0);
  for (int i = 0; i < Size; i++) {
    SumStates[0] = +X_inf[i];
    SumStates[1] = +X_sup[i];
  }
  return SumStates;
}

int main() {
  // Model to use
  RIDMME_INTERVAL_AUTODIFF(*ModelName)
  (INTERVAL_VECTOR &, INTERVAL_VECTOR &, int const) = ModelExample;
  // file names
  string FileName = "DengueDiscretedModel_exp_1";

  // Intervals for model states
  Array X_inf{332, 9.5};
  Array X_sup{367, 10.5};

  Array P_inf{0.3, 0.2};
  Array P_sup{0.3, 0.2};

  // Array P_inf{0.76, 0.297, 0.297};
  // Array P_sup{0.84, 0.3030, 0.297};

  int t_f = 5;

  ArrayInt N = {8, 8, 8, 8};

  // Series of simulations

  // Montecarlo simulations
  int Samples = 1000;
  Matrix Monte_inf(t_f, Array(X_inf.size())),
      Monte_sup(t_f, Array(X_inf.size()));
  MontecarloSets(ModelName, X_inf, X_sup, P_inf, P_sup, Samples, t_f, Monte_inf,
                 Monte_sup);
  saveMatricesToCSV(Monte_inf, Monte_sup,
                    "Results/" + FileName + "_Montecarlo.csv");

  // RIDMME simulation
  Matrix RIDMME_inf(t_f, Array(X_inf.size())),
      RIDMME_sup(t_f, Array(X_inf.size()));
  Simulation(ModelName, X_inf, X_sup, P_inf, P_sup, N, t_f, RIDMME_inf,
             RIDMME_sup);
  saveMatricesToCSV(RIDMME_inf, RIDMME_sup,
                    "Results/" + FileName + "_RIDMME.csv");

  // Classical Mixed-monotonic simulation
  Matrix MixMon_inf(t_f, Array(X_inf.size())),
      MixMon_sup(t_f, Array(X_inf.size()));
  SimulationMixMon(ModelName, X_inf, X_sup, P_inf, P_sup, t_f, MixMon_inf,
                   MixMon_sup);
  saveMatricesToCSV(MixMon_inf, MixMon_sup,
                    "Results/" + FileName + "_MixMon.csv");

  // //Estimation of error using a user-defined function
  Array Constraint = SumModel(X_inf, X_sup);
  Matrix ConstraintLimits =
      merge2Vectors(Array(t_f, Constraint[0]), Array(t_f, Constraint[1]));

  MatrixList Low(3), Up(3);
  Low[0] = RIDMME_inf;
  Low[1] = MixMon_inf;
  Low[2] = Monte_inf;

  Up[0] = RIDMME_sup;
  Up[1] = MixMon_sup;
  Up[2] = Monte_sup;

  MultipleError(SumModel, Low, Up, ConstraintLimits, 1, FileName);
  cout << endl;

  return 0;
}