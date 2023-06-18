#include "Models.h"
#include "RIDMME_Interval_Autodiff.h"

RIDMME_INTERVAL_AUTODIFF ModelFunction1(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    // RIDMME_INTERVAL_AUTODIFF x(v);
    RIDMME_INTERVAL_AUTODIFF p(P);
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
    //Parameter definition
    double beta, gamma, mu;
    beta = 0.1;
    gamma = 0.3;
    mu = 0.01;
    //Model functions
    RIDMME_INTERVAL_AUTODIFF v1 = (1-mu)*x(1)*(1-x(4)*(x(2)/(x(1)+x(2)+x(3)))) + mu*(x(1)+x(2)+x(3));                
    RIDMME_INTERVAL_AUTODIFF v2 = (1-mu)*(x(1)*x(4)*(x(2)/(x(1)+x(2)+x(3)))+x(2)*(1-gamma));
    RIDMME_INTERVAL_AUTODIFF v3 = (1-mu)*(x(3) + x(2)*gamma);
    //Save results
    vector<RIDMME_INTERVAL_AUTODIFF> w(3);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;
    //Result of the ith state
    RIDMME_INTERVAL_AUTODIFF result = w[i];
    return result;
}

RIDMME_INTERVAL_AUTODIFF ModelFunction2(INTERVAL_VECTOR & v, int const i){
    RIDMME_INTERVAL_AUTODIFF x(v);
    //Parameter definition
    double beta = 0.1, gamma = 0.3, mu = 0.01, z = 1;
    //Model functions
    RIDMME_INTERVAL_AUTODIFF G1 =  1-Power(1-Power((x(2)/(x(1)+x(2)+x(3))),(x(2)/(x(1)+x(2)))),z);
    // cout << Inf(G1.grd()) << endl;
    // RIDMME_INTERVAL_AUTODIFF G = 1-Power(x(1)/(x(2)+x(3)+x(1)), z);
    RIDMME_INTERVAL_AUTODIFF v1 = (1-mu)*x(1)*(1-beta*G1) + mu*(x(1)+x(2)+x(3));                
    RIDMME_INTERVAL_AUTODIFF v2 = (1-mu)*(x(1)*beta*G1+x(2)*(1-gamma));
    RIDMME_INTERVAL_AUTODIFF v3 = (1-mu)*(x(3) + x(2)*gamma);
    //Save results
    vector<RIDMME_INTERVAL_AUTODIFF> w(3);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;
    //Result of the ith state
    RIDMME_INTERVAL_AUTODIFF result = w[i];
    return result;
}

RIDMME_INTERVAL_AUTODIFF ModelFunction5(INTERVAL_VECTOR & v, int const i){
    RIDMME_INTERVAL_AUTODIFF x(v);
    //Parameter definition
    double beta, gamma, mu;
    beta = 0.1;
    gamma = 0.3;
    mu = 0.01;
    //Model functions
    RIDMME_INTERVAL_AUTODIFF G = 1-1/(x(2)/(x(1)+x(2)+x(3))+1);
    RIDMME_INTERVAL_AUTODIFF v1 = (1-mu)*x(1)*(1-beta*G) + mu*(x(1)+x(2)+x(3));                
    RIDMME_INTERVAL_AUTODIFF v2 = (1-mu)*(x(1)*beta*G+x(2)*(1-gamma));
    RIDMME_INTERVAL_AUTODIFF v3 = (1-mu)*(x(3) + x(2)*gamma);
    //Save results
    vector<RIDMME_INTERVAL_AUTODIFF> w(3);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;
    //Result of the ith state
    RIDMME_INTERVAL_AUTODIFF result = w[i];
    return result;
}

RIDMME_INTERVAL_AUTODIFF LogisticFunction(INTERVAL_VECTOR & v, int const i){
    RIDMME_INTERVAL_AUTODIFF x(v);

    double r = 0.35, k = 30;

    RIDMME_INTERVAL_AUTODIFF v1 =   x(1)+r*x(1)*(1-x(1)/k);
    vector<RIDMME_INTERVAL_AUTODIFF> w(1);
    w[0] = v1;
    RIDMME_INTERVAL_AUTODIFF result = w[i];

    return result;
}
