
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
#include <bits/stdc++.h>

// Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)
// Note: config/systemFolderConfi/BiasRnd.h comment line 42 to compute Nan-inf.

int main(int argc, char* argv[]) {  
  // Model to use
  int ModelNumber = std::stoi(argv[1]);
  //int ModelNumber = 4;
  cout << "Model selected: "+  to_string(ModelNumber) << endl;
  RIDMME_INTERVAL_AUTODIFF(*ModelName) (INTERVAL_VECTOR &, INTERVAL_VECTOR &, int const, int const);
  Array (*ErrorExpression)(Array X_inf, Array X_sup);
  // time coutner
  time_t start, end;
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
  time(&start);
  Matrix MixMon_inf(t_f, Array(X_inf.size())),
      MixMon_sup(t_f, Array(X_inf.size()));
  SimulationMixMon(ModelName, X_inf, X_sup, P_inf, P_sup, t_f, MixMon_inf,
                   MixMon_sup);
  time(&end);
  saveMatricesToCSV(MixMon_inf, MixMon_sup,
                    "../../Example/Results/" + FileName + "_MixMon.csv");
    // Calculating total time taken by the program.
    double time_taken = double(end - start);
    cout << "Time taken by classic MixMon is : " << fixed
        << time_taken << setprecision(5);
    cout << " sec " << endl;    

  // RIDMME simulation
  time(&start);
  Matrix RIDMME_inf(t_f, Array(X_inf.size())),
      RIDMME_sup(t_f, Array(X_inf.size()));
  Simulation(ModelName, X_inf, X_sup, P_inf, P_sup, N, t_f, RIDMME_inf,
             RIDMME_sup);
  time(&end);
  saveMatricesToCSV(RIDMME_inf, RIDMME_sup,
                    "../../Example/Results/" + FileName + "_RIDMME.csv");
    
    // Calculating total time taken by the program.
    time_taken = double(end - start);
    cout << "Time taken by RIDMME is : " << fixed
        << time_taken << setprecision(5);
    cout << " sec " << endl;      

  //Estimation of error using a user-defined function
  Array Constraint = ErrorExpression(X_inf, X_sup);
  Matrix ConstraintLimits =
      merge2Vectors(Array(t_f, Constraint[0]), Array(t_f, Constraint[1]));

  MatrixList Low(3), Up(3);
  Low[0] = RIDMME_inf; Low[1] = MixMon_inf; Low[2] = Monte_inf;

  Up[0] = RIDMME_sup; Up[1] = MixMon_sup; Up[2] = Monte_sup;

  MultipleError(ErrorExpression, Low, Up, ConstraintLimits, 1, "../../Example/Results/Error" + FileName);
  cout << endl;

  return 0.0;
}