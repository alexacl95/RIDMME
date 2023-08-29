#ifndef MO_H
#define MO_H
#include "AutomaticDifferentiation.h"
#include "Functions.h"
#include "MatrixRIDMME.h"

RIDMME_INTERVAL_AUTODIFF ModelFunction1(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i);
RIDMME_INTERVAL_AUTODIFF ModelFunction2(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i);
RIDMME_INTERVAL_AUTODIFF ModelFunction3(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i);
RIDMME_INTERVAL_AUTODIFF ModelFunction4(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i);
RIDMME_INTERVAL_AUTODIFF ModelExample(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                      int const i);
RIDMME_INTERVAL_AUTODIFF ModelDengue(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                     int const i);
RIDMME_INTERVAL_AUTODIFF ModelCOVID(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                    int const i);
RIDMME_INTERVAL_AUTODIFF ModelFunctionSin(INTERVAL_VECTOR &v, INTERVAL_VECTOR &P,
                                        int const i);                                    
RIDMME_INTERVAL_AUTODIFF DengueDiscretedModelPLOS(INTERVAL_VECTOR &v,
                                              INTERVAL_VECTOR &P, int const i);
Array SumModel(Array X_inf, Array X_sup);
Array SumModelDengue(Array X_inf, Array X_sup);
void SetModelFactors(int ModelNumber,  RIDMME_INTERVAL_AUTODIFF(**ModelName)
  (INTERVAL_VECTOR &, INTERVAL_VECTOR &, int const), Array& X_inf, Array& X_sup,
  Array& P_inf, Array& P_sup, int& t_f, ArrayInt& N, Array (**ErrorExpression)(Array X_inf, Array X_sup));
#endif