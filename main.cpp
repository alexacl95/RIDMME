
#include <stdio.h>
#include <iostream>
#include <vector>
#include "Interval.h"
#include "Functions.h"
#include "IntervalVector.h"
#include "AutomaticDifferentiation.h"
#include "Matrix.h"
#include "Refinement.h"
#include "NumericalMethod.h"
#include "Models.h"
#include "Montecarlo.h"

//Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)
//Note: config/systemFolderConfi/BiasRnd.h comment line 42 to compute Nan-inf.

int main(){
    // Intervals for model states
    array X_inf{230, 50, 20};
    array X_sup{300, 60, 40};

    // array P_inf{};
    // array P_sup{};
    
    // array P_inf{0.1};
    // array P_sup{0.3}; 

    // array P_inf{0.1, 0.3};
    // array P_sup{0.3, 0.6}; 

    array P_inf{0.1, 0.01};
    array P_sup{0.3, 0.03};

    int t_f = 25;
         
    arrayInt N = {15,15,15,2,2};

    Simulation(ModelFunction1, X_inf, X_sup, P_inf, P_sup, N, t_f);    

    int Samples = 1000;
    MontecarloSets(ModelFunction1, t_f, X_inf, X_sup, P_inf, P_sup, Samples);

    return 0;
}