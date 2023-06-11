#include "AutomaticDifferentiation.h"

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
    //Save results
    vector<INTERVAL_AUTODIFF> w(3);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;
    //Result of the ith state
    INTERVAL_AUTODIFF result = w[i];

    return result;
}

void AutoDiffModel(array X_inf, array X_sup, matrix& GrdInf, matrix& GrdSup){   
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

void MultipleAutoDiffModel(matrix X_sub, matrixInt allPerms, arrayInt N, matrix Jx_inf, matrix Jx_sup){

    int NumPers = allPerms.size();
    int n_x = allPerms[0].size();

    array X_inf(n_x);
    array X_sup(n_x);

    matrix GrdInf, GrdSup;

    N = getCumulativeSum(N);

    for (int i = 0; i < NumPers; i++){
        arrayInt row = allPerms[i];
        for(int j = 0; j < n_x; j++){
            int idx = row[j];
            X_inf[j] = X_sub[N[j]+idx-1][0];
            X_sup[j] = X_sub[N[j]+idx-1][1];
        }
        AutoDiffModel(X_inf, X_sup, GrdInf, GrdSup);
    }
}