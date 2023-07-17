#include "Montecarlo.h"

double RandN(double low, double up){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> distribution(low, up);

    double result;
    result = distribution(gen);

    return result;
}

void MontecarloSets(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const),
                                   const array& x_low, const array& x_up,
                                   const array& p_low, const array& p_up, int sampleN, int& time,
                                   matrix& Result_inf, matrix& Result_sup){
    cout << "Stariting Montecarlo" << endl;
    int n_x = x_low.size();
    int n_p = p_low.size(); 

    matrixList Montecarlo_sim(sampleN, matrix(time, array(n_x,0)));
    array x(n_x,0);
    matrix p(sampleN, array(n_p,0));
    
    for(int i = 0; i < sampleN ; i++){
        for(int j = 0; j< n_x; j++){
            x[j] = RandN(x_low[j], x_up[j]);
        }
        Montecarlo_sim[i][0] = x;
        for(int j = 0; j< n_p; j++){
            p[i][j] = RandN(p_low[j], p_up[j]);
        }
    }
    
    array x_step(n_x,0); 

    for(int i = 0; i < time-1; i++){
        // #pragma omp parallel for
        for(int j = 0;  j < sampleN; j++){
            x_step = Montecarlo_sim[j][i];
            ToDoubleEvaluationModel(fnc, i, x_step, p[j]);
            Montecarlo_sim[j][i+1] = x_step;
        }            
    }

    Result_inf = findMinMaxMatrix(Montecarlo_sim, 0);
    Result_sup = findMinMaxMatrix(Montecarlo_sim, 1);

    saveMatricesToCSV(Result_inf, Result_sup, "ResultsMontecarlo.csv");

    cout << "Finished!" << endl;
}