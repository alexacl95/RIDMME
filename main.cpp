#include <stdio.h>
#include <iostream>
#include <vector>
#include "Interval.h"
#include "Functions.h"
#include "IntervalVector.h"
#include "AutomaticDifferentiation.h"
#include "Matrix.h"
#include "Refinement.h"

//Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)

int main(){
    // Intervals for model states
    array X_inf{230, 50, 20};
    array X_sup{300, 60, 40};

    matrix X;

    MergeTwoArray(X_inf, X_sup, X);
    arrayInt N = {5,5,5};

    matrix X_sub = NonUniformSubdivision(X, N);
    matrixInt allPerms = generateAllPermutations(N);

    matrix Jx_inf(3,array(3)),Jx_sup(3,array(3));
    MultipleAutoDiffModel(X_sub, allPerms, N, Jx_inf,Jx_sup);

    PrintMatrix(Jx_inf);
    PrintMatrix(Jx_sup);

    return 0;
}