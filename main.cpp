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

//Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)

int main(){
    // Intervals for model states
    array X_inf{230, 50, 20};
    array X_sup{300, 60, 40};

    array P_inf{0.1, 0.3, 0.1};
    array P_sup{0.3, 0.4, 0.3}; 
    
    int t_f = 4;
         
    arrayInt N = {10,10,10};
    Simulation(ModelFunction1, X_inf, X_sup, P_inf, P_sup, N, t_f);
    
    return 0;
}