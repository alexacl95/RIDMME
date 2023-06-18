#include "Refinement.h"


matrix NonUniformSubdivision(array X_inf, array X_sup, array P_inf, array P_sup, arrayInt& N){
    matrix X_sub;
    matrix X;
    MergeTwoArray(X_inf, X_sup, X);
    int n_x = X.size();

    for (int i = 0; i < n_x; i++) {
        vector<double> Xi = X[i];
        double w_sub = (Xi[1] - Xi[0]) / N[i];

        for (int j = 0; j < N[i]; j++) {
            double lower_bound = Xi[0] + j * w_sub;
            double upper_bound = Xi[0] + (j + 1) * w_sub;
            X_sub.push_back({ lower_bound, upper_bound });
        }
    }
    return X_sub;
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

