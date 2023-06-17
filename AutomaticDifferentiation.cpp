#include "AutomaticDifferentiation.h"

void AutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), array X_inf, array X_sup, matrix& GrdInf, matrix& GrdSup){   
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
        Grd = GradientValue (fnc(x, i));
        GrdInf[i].resize(n_x), GrdSup[i].resize(n_x); 
        for(int j = 0; j < n_x; j++){
            GrdInf[i][j] = Inf(Grd(j+1));
            GrdSup[i][j] = Sup(Grd(j+1));
        }
    }
}

void MultipleAutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), matrix X_sub, matrixInt allPerms, arrayInt N, matrix& Jx_inf, matrix& Jx_sup){

    int NumPers = allPerms.size();
    int n_x = allPerms[0].size();

    array X_inf(n_x);
    array X_sup(n_x);

    matrix GrdInf_temp, GrdSup_temp;
    matrix GrdInf(n_x, array(n_x, numeric_limits<double>::infinity()));
    matrix GrdSup(n_x, array(n_x, -numeric_limits<double>::infinity()));

    N = getCumulativeSum(N);

    for (int i = 0; i < NumPers; i++){
        GrdInf_temp = GrdInf;
        GrdSup_temp = GrdSup;
        arrayInt row = allPerms[i];
        for(int j = 0; j < n_x; j++){
            int idx = row[j];
            X_inf[j] = X_sub[N[j]+idx-1][0];
            X_sup[j] = X_sub[N[j]+idx-1][1];
        }
        AutoDiffModel(fnc, X_inf, X_sup, GrdInf, GrdSup);
        GrdInf = combineMatrices(GrdInf_temp, GrdInf, 0);
        GrdSup = combineMatrices(GrdSup_temp, GrdSup, 1);        
    }    
    Jx_inf = GrdInf;
    Jx_sup = GrdSup;
}

void ToDoubleEvaluationModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), double t, array& x, array p){
    int n_x = x.size();
    INTERVAL_VECTOR X(n_x);
    for (int i = 0; i < n_x; i++){
        X(i+1) = INTERVAL(x[i], x[i]);
    } 
    array Sol(n_x);
    for (int i = 0; i < n_x; i++){
        Sol[i] = Inf(FunctionValue(fnc(X, i)));
    }    
    x = Sol;
}