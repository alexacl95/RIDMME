#include "Models.h"
#include "RIDMME_Interval_Autodiff.h"

RIDMME_INTERVAL_AUTODIFF ModelFunction1(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i) {
  // Martcheva-based model for a SIR with mortality model
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Parameter definition
  RIDMME_INTERVAL_AUTODIFF beta = x(4);
//   RIDMME_INTERVAL_AUTODIFF gamma = x(5);
  double gamma = 0.3;
  // Prob
  RIDMME_INTERVAL_AUTODIFF G = beta * (x(2) / (x(1) + x(2) + x(3)));
  // Model functions
  RIDMME_INTERVAL_AUTODIFF v1 = x(1) * (1 - G);
  RIDMME_INTERVAL_AUTODIFF v2 = x(1) * G + x(2) * (1 - gamma);
  RIDMME_INTERVAL_AUTODIFF v3 = x(3) + x(2) * gamma;
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(3);
  w[0] = v1;
  w[1] = v2;
  w[2] = v3;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

RIDMME_INTERVAL_AUTODIFF ModelFunction2(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i) {
  // Martcheva-based model for a SEIR model
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Parameter definition
  RIDMME_INTERVAL_AUTODIFF beta = x(5);  // x(5)
//   RIDMME_INTERVAL_AUTODIFF gamma = x(6); // x(6)
  double gamma = 0.3;
  double theta = 0.3; // x(7)
  // Probability
  RIDMME_INTERVAL_AUTODIFF G = beta * (x(2) / (x(1) + x(2) + x(3) + x(4)));
  // Model functions
  RIDMME_INTERVAL_AUTODIFF v1 = (1 - G) * x(1);
  RIDMME_INTERVAL_AUTODIFF v2 = G * x(1) + (1 - theta) * x(2);
  RIDMME_INTERVAL_AUTODIFF v3 = theta * x(2) + (1 - gamma) * x(3);
  RIDMME_INTERVAL_AUTODIFF v4 = x(4) + gamma * x(3);
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(4);
  w[0] = v1;
  w[1] = v2;
  w[2] = v3;
  w[3] = v4;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

RIDMME_INTERVAL_AUTODIFF ModelFunction3(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i) {
  // Martcheva-based model for a SARS model
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Parameter definition
  double beta, gamma, theta, alpha;
  beta = 0.8;
  gamma = 0.3;
  theta = 0.3;
  alpha = 0.1;
  // Probability
  RIDMME_INTERVAL_AUTODIFF G = beta * (x(2) / (x(1) + x(2) + x(3) + x(4)));
  // Model functions
  RIDMME_INTERVAL_AUTODIFF v1 = x(1) * (1 - G);
  RIDMME_INTERVAL_AUTODIFF v2 = x(1) * G + (1 - theta) * (1 - alpha) * x(2);
  RIDMME_INTERVAL_AUTODIFF v3 = theta * x(2) + (1 - gamma) * (1 - alpha) * x(3);
  RIDMME_INTERVAL_AUTODIFF v4 = (1 - theta) * alpha * x(2) + (1 - theta) * x(4);
  RIDMME_INTERVAL_AUTODIFF v5 = theta * x(5) + (1 - gamma) * alpha * x(3) + (1 - gamma) * x(5);
  RIDMME_INTERVAL_AUTODIFF v6 = x(6) + gamma * (x(5) + x(3));
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(6);
  w[0] = v1;
  w[1] = v2;
  w[2] = v3;
  w[3] = v4;
  w[4] = v5;
  w[5] = v6;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

RIDMME_INTERVAL_AUTODIFF ModelFunction4(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i) {
  // Martcheva-based model for a SEIR model
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Parameter definition
  double gamma_h, gamma_m, mu_m, A_m, beta_hm, mu_h;
  beta_hm = 0.3;
//   RIDMME_INTERVAL_AUTODIFF beta_hm = x(5);
  RIDMME_INTERVAL_AUTODIFF beta_mh = x(5);
//   beta_mh = 0.2;
  gamma_h = 0.3;
  gamma_m = 0.01;
  mu_h = 0.001;
  mu_m = 0.35;
  A_m = 100;
//   A_h = 0;
  // Probability
  RIDMME_INTERVAL_AUTODIFF G_mh = beta_mh * (x(4) / (x(3) + x(4)));
  RIDMME_INTERVAL_AUTODIFF G_hm = beta_hm * (x(2) / (x(1) + x(2)));
  // Model functions
  RIDMME_INTERVAL_AUTODIFF v1 = (1-mu_h)*(1 - G_mh) * x(1) + (1-mu_h)*gamma_h * x(2);
  RIDMME_INTERVAL_AUTODIFF v2 = (1-mu_h)*G_mh * x(1) +  (1-mu_h)*(1 - gamma_h) * x(2);
  RIDMME_INTERVAL_AUTODIFF v3 =
      A_m + (1 - mu_m) * (1 - G_hm) * x(3) + (1 - mu_m) * gamma_m * x(4);
  RIDMME_INTERVAL_AUTODIFF v4 =
      (1 - mu_m) * G_hm * x(3) + (1 - mu_m) * (1 - gamma_h) * x(4);
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(4);
  w[0] = v1;
  w[1] = v2;
  w[2] = v3;
  w[3] = v4;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

RIDMME_INTERVAL_AUTODIFF ModelExample(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                      int const i) {
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Model functions
  RIDMME_INTERVAL_AUTODIFF v1 =
      x(1) * (1 - x(4)) * (1 - x(3) * x(2) / (x(1) + x(2))) +
      x(4) * (x(1) + x(2));
  RIDMME_INTERVAL_AUTODIFF v2 =
      x(1) * (1 - x(4)) * (x(3) * x(2) / (x(1) + x(2))) + (1 - x(4)) * x(2);
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(2);
  w[0] = v1;
  w[1] = v2;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

RIDMME_INTERVAL_AUTODIFF ModelDengue(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                     int const i) {
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Parameters
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
  // Model functions
  RIDMME_INTERVAL_AUTODIFF aleph_H = 1 + nu_h * x(2) / (x(1) + x(2));
  RIDMME_INTERVAL_AUTODIFF aleph_M = 1 + nu_m * x(5) / (x(4) + x(5));
  RIDMME_INTERVAL_AUTODIFF red_m =
      z * x(5) * Power((x(1) / (x(1) + x(2) + x(3))), aleph_H);
  RIDMME_INTERVAL_AUTODIFF phi = beta_mh * (1 - Power(1 - 1 / x(1), red_m));
  RIDMME_INTERVAL_AUTODIFF psi =
      beta_hm * (1 - Power(1 - Power(x(2) / (x(1) + x(2) + x(3)), aleph_M), z));
  RIDMME_INTERVAL_AUTODIFF A = r * (x(4) + x(5)) * Exp(1 - (x(4) + x(5)) / C);
  // Human population
  RIDMME_INTERVAL_AUTODIFF v1 =
      x(1) * (1 - phi) * (1 - mu_h) + (x(1) + x(2) + x(3)) * mu_h;
  RIDMME_INTERVAL_AUTODIFF v2 =
      x(1) * phi * (1 - mu_h) + x(2) * (1 - gamma) * (1 - mu_h);
  RIDMME_INTERVAL_AUTODIFF v3 = x(3) * (1 - mu_h) + x(2) * gamma * (1 - mu_h);
  // Mosquitoes population
  RIDMME_INTERVAL_AUTODIFF v4 = x(4) * (1 - psi) * (1 - mu_m) + A;
  RIDMME_INTERVAL_AUTODIFF v5 = x(4) * psi * (1 - mu_m) + x(5) * (1 - mu_m);
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(5);
  w[0] = v1;
  w[1] = v2;
  w[2] = v3;
  w[3] = v4;
  w[4] = v5;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

RIDMME_INTERVAL_AUTODIFF ModelCOVID(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                    int const i) {
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Parameter definition
  double lambda = 0.022;
  double alpha = 0.047;
  double mu = 0;
  double theta = 1;
  double gamma = 0.04;
  double nu = 19;
  double beta = 0.99;
  double z = 2;
  // Model functions
  RIDMME_INTERVAL_AUTODIFF aleph =
      1 + nu * (x(3) + x(4) + x(5)) / (x(1) + x(3) + x(4) + x(5));
  RIDMME_INTERVAL_AUTODIFF G1 =
      beta *
      (1 - Power(1 - Power(x(3) / (x(1) + x(2) + x(3) + x(4) + x(5) + x(6)), aleph), z));
  // Model functions
  RIDMME_INTERVAL_AUTODIFF v1 =
      (1 - lambda) * (1 - G1) * x(1) + alpha * x(2) + mu * x(6);
  RIDMME_INTERVAL_AUTODIFF v2 = lambda * (1 - G1) * x(1) + (1 - alpha) * x(2);
  RIDMME_INTERVAL_AUTODIFF v3 =
      G1 * x(1) + (1 - lambda) * (1 - theta) * (1 - gamma) * x(3) +
      alpha * (1 - theta) * (1 - gamma) * x(4);
  RIDMME_INTERVAL_AUTODIFF v4 = lambda * (1 - gamma) * (1 - theta) * x(3) +
                                (1 - gamma) * (1 - theta) * (1 - alpha) * x(4);
  RIDMME_INTERVAL_AUTODIFF v5 = theta * (1 - gamma) * x(3) +
                                theta * (1 - gamma) * x(4) + x(5) * (1 - gamma);
  RIDMME_INTERVAL_AUTODIFF v6 = x(6) * (1 - mu) + gamma * (x(3) + x(4) + x(5));
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(6);
  w[0] = v1;
  w[1] = v2;
  w[2] = v3;
  w[3] = v4;
  w[4] = v5;
  w[5] = v6;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

RIDMME_INTERVAL_AUTODIFF DengueDiscretedModelPLOS(INTERVAL_VECTOR &v,
                                              INTERVAL_VECTOR &P, int const i) {
  RIDMME_INTERVAL_AUTODIFF x = MergeTwoIntervalArr(v, P);
  // Step
  double h = 0.1;
  // Parameter definition
  double mu_h = 0;
  double gamma = 1;
  double beta_hm = 0.87;
  double mu_m = 0.35;
  double A_m = 100;
  double tetha_m = 0.3;
  double tetha_h = 0.3;
  // Model functions
  RIDMME_INTERVAL_AUTODIFF M = x(1) + x(2) + x(3);
  RIDMME_INTERVAL_AUTODIFF H = x(4) + x(5) + x(6) + x(7);
  // Model functions
  RIDMME_INTERVAL_AUTODIFF v1 =
      x(1) + h * (A_m - beta_hm * x(1) * x(6) / H - mu_m * x(1));
  RIDMME_INTERVAL_AUTODIFF v2 =
      x(2) + h * (beta_hm * x(1) * x(6) / H - (tetha_m + mu_m) * x(2));
  RIDMME_INTERVAL_AUTODIFF v3 = x(3) + h * (tetha_m * x(2) - mu_m * x(3));
  RIDMME_INTERVAL_AUTODIFF v4 =
      x(4) + h * (mu_h * (x(5) + x(6) + x(7)) - beta_hm * x(4) * x(3) / M);
  RIDMME_INTERVAL_AUTODIFF v5 =
      x(5) + h * (beta_hm * x(4) * x(3) / M + (tetha_h + mu_h) * x(5));
  RIDMME_INTERVAL_AUTODIFF v6 =
      x(6) + h * (tetha_h * x(5) - (gamma + mu_h) * x(6));
  RIDMME_INTERVAL_AUTODIFF v7 = x(7) + h * (gamma * x(6) - mu_h * x(7));
  // Save results
  vector<RIDMME_INTERVAL_AUTODIFF> w(7);
  w[0] = v1;
  w[1] = v2;
  w[2] = v3;
  w[3] = v4;
  w[4] = v5;
  w[5] = v6;
  w[6] = v7;
  // Result of the ith state
  RIDMME_INTERVAL_AUTODIFF result = w[i];
  return result;
}

// Metric functions for population models
Array SumModel(Array X_inf, Array X_sup) {
  int Size = X_inf.size();
  Array SumStates(2, 0);
  for (int i = 0; i < Size; i++) {
    SumStates[0] = X_inf[i] + SumStates[0];
    SumStates[1] = X_sup[i] + SumStates[1];
  }
  return SumStates;
}

Array SumModelDengue(Array X_inf, Array X_sup) {
  int Size = 3;
  Array SumStates(2, 0);
  for (int i = 0; i < Size; i++) {
    SumStates[0] = +X_inf[i];
    SumStates[1] = +X_sup[i];
  }
  return SumStates;
}

// Sett initial conditions and parameters

void SetModelFactors(int ModelNumber,  RIDMME_INTERVAL_AUTODIFF(**ModelName)
  (INTERVAL_VECTOR &, INTERVAL_VECTOR &, int const), Array& X_inf, Array& X_sup,
  Array& P_inf, Array& P_sup, int& t_f, ArrayInt& N, Array (**ErrorExpression)(Array X_inf, Array X_sup)) {
    switch (ModelNumber){
    case 1:
        *ModelName = ModelFunction1;
        X_inf = {332, 9.5, 9.5}; 
        X_sup = {367, 10.5, 10.5};
        P_inf = {0.76};
        P_sup = {0.84};
        N = {15, 15, 15, 10};
        t_f = 30;
        *ErrorExpression = SumModel;
        break;
    case 2: 
        *ModelName = ModelFunction2;
        X_inf = {332, 9.5, 9.5, 9.5};
        X_sup = {367, 10.5, 10.5, 10.5};
        P_inf = {0.76};
        P_sup = {0.84};
        N = {15,15,15,15,10};
        t_f = 35;
        *ErrorExpression = SumModel;
        break;
    case 3:        
        *ModelName = ModelFunction3;
        X_inf = {332, 9.5, 9.5, 9.5, 9.5, 0};
        X_sup = {367, 10.5, 10.5, 10.5, 10.5, 0};
        P_inf = {};
        P_sup = {};
        N = {10,10,10,10,10,10};
        t_f = 35;
        *ErrorExpression = SumModel;
        break;
    case 4: 
        *ModelName = ModelFunction4;
        X_inf = {332, 9.5, 9.5, 9.5};
        X_sup = {367, 10.5, 10.5, 10.5};
        P_inf = {0.15};
        P_sup = {0.3};
        N = {15,15,15,15,10};
        t_f = 35;
        *ErrorExpression = SumModelDengue;
        break;  
    case 5:     
        *ModelName = ModelDengue;
        X_inf = {427500, 7.6, 4.75, 26628, 3.8};
        X_sup = {472500, 8.4, 5.25, 29431, 4.2};
        P_inf = {};
        P_sup = {};
        N = {15,15,15,15,15};
        t_f = 50;
        *ErrorExpression = SumModelDengue;  
        break; 
    case 6:    
        *ModelName = ModelCOVID;
        X_inf = {78850000, 0, 18, 0, 108, 15};
        X_sup = {87150000, 0, 20, 0, 120, 17};
        P_inf = {};
        P_sup = {};
        N = {10,10,15,5,10,5};
        t_f = 40;
        *ErrorExpression = SumModel;  
        break;  
    case 7: 
        *ModelName = DengueDiscretedModelPLOS;
        X_inf = {};
        X_sup = {};
        P_inf = {};
        P_sup = {};
        t_f = 30;
        *ErrorExpression = SumModelDengue;  
        break;   
    default:
        break;
    }

}