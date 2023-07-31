#include "Models.h"
#include "RIDMME_Interval_Autodiff.h"

RIDMME_INTERVAL_AUTODIFF ModelFunction1(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
    //Parameter definition
    double beta, gamma, mu;
    beta = 0.1;
    gamma = 0.3;
    mu = 0.01;
    //Model functions
    RIDMME_INTERVAL_AUTODIFF v1 = (1-x(5))*x(1)*(1-x(4)*(x(2)/(x(1)+x(2)+x(3)))) + x(5)*(x(1)+x(2)+x(3));                
    RIDMME_INTERVAL_AUTODIFF v2 = (1-x(5))*(x(4)*x(1)*(x(2)/(x(1)+x(2)+x(3)))+x(2)*(1-gamma));
    RIDMME_INTERVAL_AUTODIFF v3 = (1-x(5))*(x(3) + x(2)*gamma);
    //Save results
    vector<RIDMME_INTERVAL_AUTODIFF> w(3);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;
    //Result of the ith state
    RIDMME_INTERVAL_AUTODIFF result = w[i];
    return result;
} 

RIDMME_INTERVAL_AUTODIFF ModelFunction2(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
    //Parameter definition
    double beta = 0.1, gamma = 0.3, mu = 0.01, z = 1;
    //Model functions
    RIDMME_INTERVAL_AUTODIFF G1 =  1-Power(1-Power((x(2)/(x(1)+x(2)+x(3))),(x(2)/(x(1)+x(2)))),z);
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

RIDMME_INTERVAL_AUTODIFF ModelFunction5(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
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

RIDMME_INTERVAL_AUTODIFF LogisticFunction(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);

    double r = 0.35, k = 30;

    RIDMME_INTERVAL_AUTODIFF v1 =   x(1)+r*x(1)*(1-x(1)/k);
    vector<RIDMME_INTERVAL_AUTODIFF> w(1);
    w[0] = v1;
    RIDMME_INTERVAL_AUTODIFF result = w[i];

    return result;
}

RIDMME_INTERVAL_AUTODIFF ModelExample(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
    //Model functions
    RIDMME_INTERVAL_AUTODIFF v1 = x(1)*(1-x(4))*(1-x(3)*x(2)/(x(1)+x(2))) + x(4)*(x(1)+x(2));                
    RIDMME_INTERVAL_AUTODIFF v2 = x(1)*(1-x(4))*(x(3)*x(2)/(x(1)+x(2))) + (1-x(4))*x(2);
    //Save results
    vector<RIDMME_INTERVAL_AUTODIFF> w(2);
    w[0] = v1;
    w[1] = v2;
    //Result of the ith state
    RIDMME_INTERVAL_AUTODIFF result = w[i];
    return result;
}

RIDMME_INTERVAL_AUTODIFF ModelDengue(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
    //Parameters
    double mu_h = 0;
    double gamma = 1;
    double beta_hm = 0.87;
    double beta_mh = 0.28;
    double mu_m = 0.35;
    double z = 8;
    double C = 325000;
    double r = 0.12;
    double nu_h = 0.0006;
    double nu_m = 0.8;
    //Model functions
    RIDMME_INTERVAL_AUTODIFF aleph_H = 1 + nu_h*x(2)/(x(1)+x(2));
    RIDMME_INTERVAL_AUTODIFF aleph_M = 1 + nu_m*x(5)/(x(4)+x(5));
    RIDMME_INTERVAL_AUTODIFF red_m = z*x(5)*Power((x(1)/(x(1)+x(2)+x(3))), aleph_H);
    RIDMME_INTERVAL_AUTODIFF phi = beta_mh*(1-Power((1-1/x(1)), red_m));
    RIDMME_INTERVAL_AUTODIFF psi = beta_hm*(1-Power(1-Power(x(2)/(x(1)+x(2)+x(3)), aleph_M), z));
    RIDMME_INTERVAL_AUTODIFF A = r *(x(4)+x(5))*Exp(1-(x(4)+x(5))/C);
    //Human population
    RIDMME_INTERVAL_AUTODIFF v1 = x(1)*(1-phi)*(1-mu_h) + (x(1)+x(2)+x(3))*mu_h;
    RIDMME_INTERVAL_AUTODIFF v2 = x(1)*phi*(1-mu_h) + x(2)*(1-gamma)*(1-mu_h);
    RIDMME_INTERVAL_AUTODIFF v3 = x(3)*(1-mu_h) + x(2)*gamma*(1-mu_h) ;
    //Mosquitoes population
    RIDMME_INTERVAL_AUTODIFF v4 = x(4)*(1-psi)*(1-mu_m) + A;
    RIDMME_INTERVAL_AUTODIFF v5 = x(4)*psi*(1-mu_m) + x(5)*(1-mu_m);
    //Save results
    vector<RIDMME_INTERVAL_AUTODIFF> w(5);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;
    w[3] = v4;
    w[4] = v5;
    //Result of the ith state
    RIDMME_INTERVAL_AUTODIFF result = w[i];
    return result;
}

RIDMME_INTERVAL_AUTODIFF ModelCOVID(INTERVAL_VECTOR & v, INTERVAL_VECTOR& P, int const i){
    RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
    //Parameter definition
    double lambda = 0.022;
    double alpha = 0.047;
    double mu = 0.0;
    double theta = 1;
    double gamma = 0.04;
    double nu = 19;
    double beta = 0.99;
    double z = 2;
    //Model functions
    RIDMME_INTERVAL_AUTODIFF aleph = 1 + nu*(x(3)+x(4)+x(5))/(x(1) + x(3)+x(4)+x(5));
    RIDMME_INTERVAL_AUTODIFF G1 = beta*(1-Power(1-Power(x(3)/(x(1)+x(2)+x(3)+x(4)+x(5)),aleph),z));
    //Model functions
    RIDMME_INTERVAL_AUTODIFF v1 = (1-lambda)*(1-G1)*x(1) + alpha*x(2) + mu*x(6);
    RIDMME_INTERVAL_AUTODIFF v2 = lambda*(1-G1)*x(1) + (1-alpha)*x(2);
    RIDMME_INTERVAL_AUTODIFF v3 = G1*x(1) + (1-lambda)*(1-theta)*(1-gamma)*x(3) + alpha*(1 - theta)*(1-gamma)*x(4);
    RIDMME_INTERVAL_AUTODIFF v4 = lambda*(1-gamma)*(1 - theta)*x(3) + (1 - gamma)*(1-theta)*(1-alpha)*x(4);
    RIDMME_INTERVAL_AUTODIFF v5 = theta*(1-gamma)*x(3) + theta*(1 - gamma)*x(4) + x(5)*(1-gamma);
    RIDMME_INTERVAL_AUTODIFF v6 = x(6)*(1-mu) + gamma*(x(3) + x(4) + x(5));
    //Save results
    vector<RIDMME_INTERVAL_AUTODIFF> w(6);
    w[0] = v1;
    w[1] = v2;
    w[2] = v3;
    w[3] = v4;
    w[4] = v5;
    w[5] = v6;
    //Result of the ith state
    RIDMME_INTERVAL_AUTODIFF result = w[i];
    return result;
} 