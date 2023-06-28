#include "NumericalMethod.h"

void MixMonotonicityMethod(
    RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const),
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

void Simulation(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const), array X_inf, array X_sup, array P_inf, array P_sup, arrayInt N, int t_f){
    //Check that user 
    CheckSimSettings(fnc, X_inf, X_sup, P_inf, P_sup, N);

    int n_x = X_inf.size();
    int n_p = P_inf.size();
    int n_f = n_x + n_p;
    arrayInt IterationRep(t_f);
    arrayInt N0 = N;
    int Step = 3;
    int Step0 = Step;

    matrix Results_low(t_f,array(n_x)), Results_up(t_f,array(n_x));
    // array P_inf(n_x), P_sup(n_x);

    Results_low[0] = X_inf;
    Results_up[0] = X_sup;

    //find if Any Element is a Nan 
    vector <bool> AnyNanElement(n_x,false);   

    for (int i = 1; i < t_f; i++){
        cout << "Iteration: " << i << endl;

        SingleSimulation(fnc, X_inf, X_sup, P_inf, P_sup, N, t_f, i, AnyNanElement);
        
        Results_low[i] = X_inf;
        Results_up[i] = X_sup;

        //Adapatative Refinement
        
        int Counter = 0;
        
        ReduceN(N, N0, 3, 0.8);

        while(any_of(AnyNanElement.begin(), AnyNanElement.end(), [](bool value) { return value; }) && Counter < 2){
            cout << "Starting refinement process" << endl;

            RefinementIteration(X_inf, X_sup, Results_low, Results_up, N, i,
                        IterationRep, Counter, Step, 
                        Step0, AnyNanElement, n_x, n_p);

            SingleSimulation(fnc, X_inf, X_sup, P_inf, P_sup, N, t_f, i, AnyNanElement);  

            Results_low[i] = X_inf;
            Results_up[i] = X_sup;            
        }

        IterationRep[i] =+ 1;
        
        //Finishing Refinement
    }
    //Save results in a csv file
    saveMatricesToCSV(Results_low, Results_up, "Results.csv");
}

int CheckSimSettings(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const),
                    array X_inf,array X_sup, array P_inf, array P_sup, arrayInt N){

    if(X_inf.size()!= X_sup.size()){
        throw invalid_argument("X_low and X_up arrays are not the same length.\nCheck the superior and inferior limits for model states");
        return -1;
    }else if (P_inf.size()!= P_sup.size()){
        throw invalid_argument("P_low and P_up arrays are not the same length.\nCheck the superior and inferior limits for model parameters");
        return -1;
    }

    if(N.size() != X_inf.size() + P_inf.size()){
        throw invalid_argument("Size of subdivision array is not the same as the size of factors (states+parameters).\nCheck: assing one element of N per each state/parameter.");
        return -1;
    }

    return 0;    
}

void SingleSimulation(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR&, INTERVAL_VECTOR&, int const),
                        array& X_inf, array& X_sup, array P_inf, array P_sup, arrayInt N, int t_f, int i, arrayBool& AnyNanElement){

    int n_x = X_inf.size();
    int n_p = P_inf.size();
    int n_f = n_x + n_p;
    array Factor_inf(n_f), Factor_sup(n_f);                        

    //Subdivision process of the states
    AppTwoArray(X_inf, P_inf, Factor_inf);
    AppTwoArray(X_sup, P_sup, Factor_sup);
    matrix X_sub = NonUniformSubdivision(Factor_inf, Factor_sup, N);
    matrixInt allPerms = generateAllPermutations(N);
    cout << "Number of subdivisions:"<< allPerms.size() <<endl;
           
    //Calculation of refined Jacobian
    matrix Jx_inf(n_x,array(n_x)), Jx_sup(n_x,array(n_x));
    matrix Jp_inf(n_x,array(n_p)),Jp_sup(n_x,array(n_p));

    MultipleAutoDiffModel(fnc, X_sub, allPerms, N, Jx_inf, Jx_sup, Jp_inf, Jp_sup, AnyNanElement);
    //mixed monotonicity numerical method;
    MixMonotonicityMethod(fnc, i,
                    X_inf, X_sup, P_inf, P_sup, 
                    Jx_inf,Jx_sup,Jp_inf,Jp_sup);
    cout<<endl;
}
