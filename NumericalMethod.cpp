#include "NumericalMethod.h"

void MixMonotonicityMethod(
    RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const),
    double t,
    array& x_low, array& x_up,
    array p_low, array p_up,
    matrix J_x_low, matrix J_x_up,
    matrix J_p_low, matrix J_p_up
){
    int n_x = x_low.size();
    int n_p = p_low.size();

    if (any_of(p_low.begin(), p_low.end(), [](double val) { return isnan(val); }) ||
        any_of(p_up.begin(), p_up.end(), [](double val) { return isnan(val); })) {
        p_low.assign(n_x, 0.0);
        p_up.assign(n_x, 0.0);
    }

    matrix xi_low(n_x, array(n_x));
    matrix xi_up(n_x, array(n_x));
    matrix comp_state(n_x, array(n_x));
    matrix pi_low(n_x, array(n_p));
    matrix pi_up(n_x, array(n_p));
    matrix comp_input(n_x, array(n_p));

    array Complement(n_x);
    array ComplementStates(n_x);
    array ComplementInputs(n_x);

    array succ_low(n_x);
    array succ_up(n_x);

    for (int i = 0; i < n_x; i++) {
        // Define the state auxiliary variables based on the state Jacobian
        for (int j = 0; j < n_x; j++) {
            if ((J_x_low[i][j] + J_x_up[i][j]) / 2.0 >= 0) {
                // Center of Jacobian interval is positive
                xi_low[i][j] = x_low[j];
                xi_up[i][j] = x_up[j];
                if (J_x_low[i][j] < 0) {
                    // Internal not fully positive
                    comp_state[i][j] = J_x_low[i][j];
                }
            } else {
                // Center of Jacobian interval is negative
                xi_low[i][j] = x_up[j];
                xi_up[i][j] = x_low[j];
                if (J_x_up[i][j] > 0) {
                    // Internal not fully negative
                    comp_state[i][j] = J_x_up[i][j];
                }
            }
        }
        // Define the input auxiliary variables based on the input Jacobian
        for (int j = 0; j < n_p; j++) {
            if ((J_p_low[i][j] + J_p_up[i][j]) / 2.0 >= 0) {
                // Center of Jacobian interval is positive
                pi_low[i][j] = p_low[j];
                pi_up[i][j] = p_up[j];
                if (J_p_low[i][j] < 0) {
                    // Internal not fully positive
                    comp_input[i][j] = J_p_low[i][j];
                }
            } else {
                // Center of Jacobian interval is negative
                pi_low[i][j] = p_up[j];
                pi_up[i][j] = p_low[j];
                if (J_p_up[i][j] > 0) {
                    // Internal not fully negative
                    comp_input[i][j] = J_p_up[i][j];
                }
            }
        }
        array xi_low_temp = xi_low[i];
        array xi_up_temp = xi_up[i];
        
        ToDoubleEvaluationModel(fnc,t, xi_low_temp, pi_low[i]);
        array xi_low_succ = xi_low_temp;
        ToDoubleEvaluationModel(fnc, t, xi_up_temp, pi_up[i]);
        array xi_up_succ = xi_up_temp;

        // Over-approximation on dimension i
        ComplementStates = matrixVectorMultiply(comp_state,susVectors(xi_low[i],xi_up[i]));
        ComplementInputs = matrixVectorMultiply(comp_input,susVectors(pi_low[i],pi_up[i]));
        Complement = susVectors(ComplementStates, ComplementInputs);
        succ_low[i] = xi_low_succ[i] - Complement[i];
        succ_up[i] = xi_up_succ[i] + Complement[i];
    }

    x_low = succ_low;
    x_up = succ_up;

    PrintArray(x_low);
    PrintArray(x_up);
}

void Simulation(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, int const), array X_inf, array X_sup, arrayInt N, int t_f){
    
    int n_x = X_inf.size();
    
    matrix Results_low(t_f,array(n_x)), Results_up(t_f,array(n_x));
    matrix Jp_inf(n_x,array(n_x)),Jp_sup(n_x,array(n_x));
    array P_inf(n_x), P_sup(n_x);

    Results_low[0] = X_inf;
    Results_up[0] = X_sup;

    for (int i = 1; i < t_f; i++){
        cout << "Iteration: " << i << endl;

        //Subdivision process of the states
        matrix X_sub = NonUniformSubdivision(X_inf, X_sup, N);
        matrixInt allPerms = generateAllPermutations(N);
        cout << "Number of subdivisions:"<< allPerms.size() <<endl;
        
        //Calculation of refined Jacobian
        matrix Jx_inf(n_x,array(n_x)), Jx_sup(n_x,array(n_x));
        MultipleAutoDiffModel(fnc, X_sub, allPerms, N, Jx_inf, Jx_sup);
       
        //mixed monotonicity numerical method
        MixMonotonicityMethod(fnc, i,
                        X_inf, X_sup, P_inf, P_sup, 
                        Jx_inf,Jx_sup,Jp_inf,Jp_sup);
        cout << endl; 
        Results_low[i] = X_inf;
        Results_up[i] = X_sup;
    }
    //Save results in a csv file
    saveMatricesToCSV(Results_low, Results_up, "Results.csv");
}