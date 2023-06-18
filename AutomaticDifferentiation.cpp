#include "AutomaticDifferentiation.h"

void AutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&,INTERVAL_VECTOR&, int const),
                    array X_inf, array X_sup, array P_inf, array P_sup, matrix& GrdInf, matrix& GrdSup, matrix& GrdInfP, matrix& GrdSupP){
    int n_x = X_inf.size();
    INTERVAL_VECTOR x(n_x);
    INTERVAL_VECTOR Grd(n_x);

    GrdInf.resize(n_x), GrdSup.resize(n_x);
    //Save inf and sup state intervals into INTERVAL
    for (int i = 0; i < n_x; i++){
        x(i+1) = INTERVAL(X_inf[i], X_sup[i]);
    }    

    int n_p = P_inf.size();
    INTERVAL_VECTOR p(n_p);
    INTERVAL_VECTOR GrdP(n_x);
    
    if(n_p!=0){
        INTERVAL_VECTOR p(n_p);
        INTERVAL_VECTOR GrdP(n_x);

        GrdInfP.resize(n_x), GrdSupP.resize(n_x);
        //Save inf and sup state intervals into INTERVAL
        for (int i = 0; i < n_p; i++){
            p(i+1) = INTERVAL(P_inf[i], P_sup[i]);
        }  
            //Calculation of inf and sup Jacobian values
        for (int i = 0; i < n_x; i++){
            GrdP = GradientValue(fnc(x, p, i));
            GrdInfP[i].resize(n_p), GrdSupP[i].resize(n_p); 
            for(int j = 0; j < n_p; j++){
                GrdInfP[i][j] = Inf(GrdP(j+1));
                GrdSupP[i][j] = Sup(GrdP(j+1));
            }
        }        
    }else{
        INTERVAL_VECTOR p(n_x);
        matrix zeros(n_x);
        GrdInfP = zeros;
        GrdSupP = zeros;
    }
    for (int i = 0; i < n_x; i++){
        Grd = GradientValue(fnc(x, p, i));
        GrdInf[i].resize(n_x), GrdSup[i].resize(n_x); 
        for(int j = 0; j < n_x; j++){
            GrdInf[i][j] = Inf(Grd(j+1));
            GrdSup[i][j] = Sup(Grd(j+1));
        }
    }
}

void MultipleAutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const), matrix X_sub, matrixInt allPerms, arrayInt N,
                        matrix& Jx_inf, matrix& Jx_sup, matrix& Jp_inf, matrix& Jp_sup){

    int NumPers = allPerms.size();
    int n_x = Jx_inf[0].size();
    int n_p = Jp_inf[0].size();

    array Factors_inf(NumPers);
    array Factors_sup(NumPers);

    array X_inf(n_x);
    array X_sup(n_x);

    array P_inf(n_p);
    array P_sup(n_p);

    matrix GrdInf_temp, GrdSup_temp;
    matrix GrdInf(n_x, array(n_x, numeric_limits<double>::infinity()));
    matrix GrdSup(n_x, array(n_x, -numeric_limits<double>::infinity()));

    matrix GrdInfP(n_x, array(n_p, numeric_limits<double>::infinity()));
    matrix GrdSupP(n_x, array(n_p, -numeric_limits<double>::infinity()));

    N = getCumulativeSum(N);

    for (int i = 0; i < NumPers; i++){
        GrdInf_temp = GrdInf;
        GrdSup_temp = GrdSup;
        arrayInt row = allPerms[i];
        for(int j = 0; j < n_x+n_p; j++){
            int idx = row[j];
            if(j<n_x){
                X_inf[j] = X_sub[N[j]+idx-1][0];
                X_sup[j] = X_sub[N[j]+idx-1][1];
            }else{
                P_inf[j-n_x] = X_sub[N[j]+idx-1][0];
                P_sup[j-n_x] = X_sub[N[j]+idx-1][1];
            }            
        }
        AutoDiffModel(fnc, X_inf, X_sup, P_inf, P_sup, GrdInf, GrdSup, GrdInfP, GrdSupP);
        GrdInf = combineMatrices(GrdInf_temp, GrdInf, 0);
        GrdSup = combineMatrices(GrdSup_temp, GrdSup, 1);        
    }    
    Jx_inf = GrdInf;
    Jx_sup = GrdSup;
}

void ToDoubleEvaluationModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const), double t, array& x, array& p){
    int n_x = x.size();
    INTERVAL_VECTOR X(n_x);
    for (int i = 0; i < n_x; i++){
        X(i+1) = INTERVAL(x[i], x[i]);
    } 
    
    int n_p = p.size();
    INTERVAL_VECTOR P(n_p);
    for (int i = 0; i < n_p; i++){
        P(i+1) = INTERVAL(p[i], p[i]);
    } 

    array Sol(n_x);
    for (int i = 0; i < n_x; i++){
        Sol[i] = Inf(FunctionValue(fnc(X,P,i)));
    }    
    x = Sol;
}