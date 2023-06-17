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
    
    int t_f = 2;
         
    arrayInt N = {1,1,1};
    Simulation(ModelFunction2, X_inf, X_sup, N, t_f);
    
    return 0;
}