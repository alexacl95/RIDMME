
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
#include <time.h>

// Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)
// Note: config/systemFolderConfi/BiasRnd.h comment line 42 to compute Nan-inf.

int main(int argc, char* argv[]) {  
  // Model to use
  int ModelNumber = std::stoi(argv[1]);
//   int ModelNumber = 2;
  cout << "Model selected: "+  to_string(ModelNumber) << endl;
  RIDMME_INTERVAL_AUTODIFF(*ModelName) (INTERVAL_VECTOR &, INTERVAL_VECTOR &, int const);
  Array (*ErrorExpression)(Array X_inf, Array X_sup);
  // file names
  string FileName = "Model_"+ to_string(ModelNumber);
  // Intervals for model states
  Array X_inf;
  Array X_sup;
  Array P_inf;
  Array P_sup;
  int t_f;
  ArrayInt N;
  SetModelFactors(ModelNumber, &ModelName, X_inf, X_sup, P_inf, P_sup, t_f, N, &ErrorExpression);

  // Series of simulations
  // Montecarlo simulations
  int Samples = 1000;
  Matrix Monte_inf(t_f, Array(X_inf.size())),
      Monte_sup(t_f, Array(X_inf.size()));
  MontecarloSets(ModelName, X_inf, X_sup, P_inf, P_sup, Samples, t_f, Monte_inf,
                 Monte_sup);
  saveMatricesToCSV(Monte_inf, Monte_sup,
                    "../../Example/Results/" + FileName + "_Montecarlo.csv");

  // Classical Mixed-monotonic simulation
  Matrix MixMon_inf(t_f, Array(X_inf.size())),
      MixMon_sup(t_f, Array(X_inf.size()));
  SimulationMixMon(ModelName, X_inf, X_sup, P_inf, P_sup, t_f, MixMon_inf,
                   MixMon_sup);
  saveMatricesToCSV(MixMon_inf, MixMon_sup,
                    "../../Example/Results/" + FileName + "_MixMon.csv");

  // RIDMME simulation
  clock_t timeSim = clock();
  Matrix RIDMME_inf(t_f, Array(X_inf.size())),
      RIDMME_sup(t_f, Array(X_inf.size()));
  Simulation(ModelName, X_inf, X_sup, P_inf, P_sup, N, t_f, RIDMME_inf,
             RIDMME_sup);
  timeSim = clock() - timeSim;
  saveMatricesToCSV(RIDMME_inf, RIDMME_sup,
                    "../../Example/Results/" + FileName + "_RIDMME.csv");
  cout << "It took " << (float(timeSim))/CLOCKS_PER_SEC << " seconds for RIDMME.\n" << endl;      

  //Estimation of error using a user-defined function
  Array Constraint = ErrorExpression(X_inf, X_sup);
  Matrix ConstraintLimits =
      merge2Vectors(Array(t_f, Constraint[0]), Array(t_f, Constraint[1]));

  MatrixList Low(3), Up(3);
  Low[0] = RIDMME_inf; Low[1] = MixMon_inf; Low[2] = Monte_inf;

  Up[0] = RIDMME_sup; Up[1] = MixMon_sup; Up[2] = Monte_sup;

  MultipleError(ErrorExpression, Low, Up, ConstraintLimits, 1, "../../Example/Results/" + FileName);
  cout << endl;

  return 0.0;
}