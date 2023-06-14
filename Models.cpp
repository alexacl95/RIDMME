#include "Models.h"

INTERVAL_AUTODIFF ModelFunction1(INTERVAL_VECTOR & v, int const i){
    INTERVAL_AUTODIFF x(v);
    //Parameter definition
    double beta, gamma, mu;
    beta = 0.1;
    gamma = 0.3;
    mu = 0.01;
    //Model functions
    INTERVAL_AUTODIFF v1 = (1-mu)*x(1)*(1-beta*(x(2)/(x(1)+x(2)+x(3)))) + mu*(x(1)+x(2)+x(3));                
    INTERVAL_AUTODIFF v2 = (1-mu)*(x(1)*beta*(x(2)/(x(1)+x(2)+x(3)))+x(2)*(1-gamma));
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

INTERVAL_AUTODIFF LogisticFunction(INTERVAL_VECTOR & v, int const i){
    INTERVAL_AUTODIFF x(v);

    double r = 0.35, k = 30;

    INTERVAL_AUTODIFF v1 =   x(1)+r*x(1)*(1-x(1)/k);
    vector<INTERVAL_AUTODIFF> w(1);
    w[0] = v1;
    INTERVAL_AUTODIFF result = w[i];

    return result;
}
