#include "Refinement.h"


matrix NonUniformSubdivision(array Factor_inf, array Factor_sup, arrayInt& N){
    matrix Factor_sub;
    matrix X;
    MergeTwoArray(Factor_inf, Factor_sup, X);
    int n_x = X.size();
    for (int i = 0; i < n_x; i++) {
        vector<double> Xi = X[i];
        double w_sub = (Xi[1] - Xi[0]) / N[i];

        for (int j = 0; j < N[i]; j++) {
            double lower_bound = Xi[0] + j * w_sub;
            double upper_bound = Xi[0] + (j + 1) * w_sub;
            Factor_sub.push_back({ lower_bound, upper_bound });
        }
    }
    return Factor_sub;
}

matrixInt generateAllPermutations(arrayInt& N) {
    int n_x = N.size();
    int totalPermutations = accumulate(N.begin(), N.end(), 1, multiplies<int>());

    matrixInt allPerms(totalPermutations, arrayInt(n_x));
    arrayInt indices(n_x);

    for (int i = 0; i < totalPermutations; i++) {
        for (int j = 0; j < n_x; j++) {
            indices[j] = (i / accumulate(N.begin() + j + 1, N.end(), 1, multiplies<int>())) % N[j] + 1;
        }
        allPerms[i] = indices;
    }
    return allPerms;
}

void ReduceN(arrayInt& N, arrayInt& N0, int TimesN0, double PercentN){
    arrayInt idx_N;
    for (size_t i = 0; i < N.size(); ++i) {
        if (N[i] >= TimesN0 * N0[i]) {
            idx_N.push_back(i);
        }
    }
    if (!idx_N.empty()) {
        for (int idx : idx_N) {
            N[idx] = ceil(PercentN * N[idx]);
        }
    cout<<"Reducing N" <<endl;
    }
}

void RefinementIteration(arrayInt& N, int& i,
                         arrayInt IterationRep, int& Step, 
                         int Step0, arrayBool AnyNanElement, int n_x, int n_p ){



    if( i != 1 & IterationRep[1] == 0){
    //Go back to first iteration and increase de N size
        i = 1;
        for (int l = 0; l <= n_x+n_p; l++){
            if(N[l] > 2*(n_x+n_p)){
                N[l] = 3*N[l];
            }else{
                N[l] = 2*(n_x+n_p);
            }
        }
        cout << "~~~~~~~~~~~~~Back to Iteration: " << i << "~~~~~~~~~~~~~" << endl;
    } else if (i > Step) {
        i -= Step;

        if(IterationRep[i] > 2){
            Step = Step+2;
            i -= Step;
            if( i < 1){
                i = 1;
            }
        }else{
            Step = Step0;
        }        

        cout << "~~~~~~~~~~~~~Back to Iteration: " << i << "~~~~~~~~~~~~~" << endl;
        for (int j = 0; j < n_x; ++j) {
            if (AnyNanElement[j]) {
                N[j] += ceil((n_x + n_p)/(Step/2));
            }
        }
        for (int j = n_x; j < n_x + n_p; ++j) {
            N[j] += ceil((n_x + n_p)/(Step/2));
        }        
    }

    cout << "  ***N = [";
        for (int val : N) {
            cout << val << " ";
        }
    cout << "]***" << endl;    
}
