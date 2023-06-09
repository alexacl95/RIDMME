#include <stdio.h>
#include "Interval.h"
#include "Functions.h"
#include "IntervalVector.h"
#include <iostream>
#include <vector>
#include "AutomaticDifferentiation.h"
#include "Matrix.h"

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
    
    return 0;
}