#include <stdio.h>
#include "Interval.h"
#include "Functions.h"
#include "IntervalVector.h"
#include <iostream>
#include <vector>
#include "AutoDiff/IntervalAutoDiff.h"
#include "AutoDiff/AutoDiff.h"

//Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)

using namespace std;

using array = vector<double>;
using matrix = vector<array>;

INTERVAL_AUTODIFF ModelFunction (INTERVAL_VECTOR & v, int const i){
    INTERVAL_AUTODIFF x(v);

    //Parameter definition
    double beta, gamma, mu;
    beta = 0.1;
    gamma = 0.3;
    mu = 0.01;

    //Model functions
    INTERVAL_AUTODIFF v1 = (1-mu)*x(1)*(1-beta*x(2)/(x(1)+x(2)+x(3))) + mu*(x(1)+x(2)+x(3));                
    INTERVAL_AUTODIFF v2 = (1-mu)*(x(1)*beta*x(2)/(x(1)+x(2)+x(3))+x(2)*(1-gamma));
    INTERVAL_AUTODIFF v3 = (1-mu)*(x(3) + x(2)*gamma);

    vector<INTERVAL_AUTODIFF> w(3);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;

    //Result of the ith state
    INTERVAL_AUTODIFF result = w[i];

    return result;
}

void AutomaticDifferentiation(array X_inf, array X_sup, matrix& GrdInf, matrix& GrdSup){   
    int n_x = X_inf.size();
    INTERVAL_VECTOR x(n_x);
    INTERVAL_VECTOR Grd(n_x);

    GrdInf.resize(n_x), GrdSup.resize(n_x);

    //Save inf and sup state intervals into INTERVAL
    for (int i = 0; i < n_x; i++){
        x(i+1) = INTERVAL(X_inf[i], X_sup[i]);
    }
    
    //Calculation of inf and sup Jacobian values
    for (int i = 0; i < n_x; i++){
        Grd = GradientValue (ModelFunction (x, i));
        GrdInf[i].resize(n_x), GrdSup[i].resize(n_x); 
        for(int j = 0; j < n_x; j++){
            GrdInf[i][j] = Inf(Grd(j+1));
            GrdSup[i][j] = Sup(Grd(j+1));
        }
    }
}

void PrintMatrix(matrix& vec){
    for (int i = 0; i < vec.size(); i++){
        for (int j = 0; j < vec[i].size(); j++){
            cout << vec[i][j]<< " ";
    }
    cout << endl;
}
}

int main()
{
    array X_inf{230, 50, 20};
    array X_sup{300, 60, 40};

    matrix GrdInf, GrdSup;
    AutomaticDifferentiation(X_inf,X_sup, GrdInf, GrdSup);

    cout << "J_x inf" << endl;
    PrintMatrix(GrdInf);
    cout << "J_x sup" << endl;
    PrintMatrix(GrdSup);
    

    
    
    
    // INTERVAL_VECTOR x1(3), x2(3);
    // x1(1) = INTERVAL(230, 300);
    // x1(2) = INTERVAL(50, 60);
    // x1(3) = INTERVAL(20, 40);

    // int idx = 2;
    // cout << "f(X) = " << FunctionValue (ModelFunction (x1,idx)) << endl;
    // cout << "inf[f'(X)] = " << Inf(GradientValue (ModelFunction (x1,idx))) << endl;
    // cout << "sup[f'(X)] = " << Sup(GradientValue (ModelFunction (x1,idx))) << endl;
    return 0;
}