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
    array X_inf{230, 50, 20};
    array X_sup{300, 60, 40};

    matrix GrdInf, GrdSup;
    AutoDiffModel(X_inf,X_sup, GrdInf, GrdSup);

    cout << "J_x inf" << endl;
    PrintMatrix(GrdInf);
    cout << "J_x sup" << endl;
    PrintMatrix(GrdSup);
    
    matrix X = {{1.0, 5.0}, {10, 20}, {1, 2}};
    vector<int> N = {2,2,4};

    matrix X_sub = NonUniformSubdivision(X, N);
    matrixInt allPerms = generateAllPermutations(N);

    PrintMatrix(X_sub);
    PrintMatrixInt(allPerms);

    return 0;
}