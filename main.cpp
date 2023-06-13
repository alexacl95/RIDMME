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

//Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)

int main(){
    // Intervals for model states
    array X_inf{230, 50, 20};
    array X_sup{300, 60, 40};

    matrix X;

    MergeTwoArray(X_inf, X_sup, X);
    arrayInt N = {30,30,30};

    matrix X_sub = NonUniformSubdivision(X, N);
    matrixInt allPerms = generateAllPermutations(N);

    matrix Jx_inf(3,array(3)),Jx_sup(3,array(3));
    MultipleAutoDiffModel(X_sub, allPerms, N, Jx_inf,Jx_sup);

    cout<< allPerms.size() <<endl;

    PrintMatrix(Jx_inf);
    PrintMatrix(Jx_sup);


    // matrix Jx_i ={{0.987316, -0.0797188, 0.0171156},
    //                 {-0.007425, 0.730125, -0.0198},
    //                 {0, 0.297, 0.99}};
    
    // matrix Jx_u = {{1.00743, -0.031745, 0.0298},
    //                 {0.0126844, 0.784884, -0.00711562},
    //                 {0, 0.297, 0.99}}; 

    //ToDoubleEvaluationModel(1.0, X_inf, X_inf);

    matrix Jp_inf(3,array(3)),Jp_sup(3,array(3));
    array P_inf(3), P_sup(3);

    PrintArray(X_inf);
    PrintArray(X_sup);
    array D = MonotonicityModel(ToDoubleEvaluationModel, 3.0,
            X_inf, X_sup, P_inf, P_sup, 
            Jx_inf,Jx_sup,Jp_inf,Jp_sup);

    cout<<endl;

    return 0;
}