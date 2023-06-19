#ifndef MO_H
#define MO_H
#include "Matrix.h"
#include "Functions.h"
#include "AutomaticDifferentiation.h"

RIDMME_INTERVAL_AUTODIFF ModelFunction1(INTERVAL_VECTOR & v, INTERVAL_VECTOR & P, int const i);
RIDMME_INTERVAL_AUTODIFF ModelFunction2(INTERVAL_VECTOR & v, INTERVAL_VECTOR & P, int const i);
RIDMME_INTERVAL_AUTODIFF ModelFunction5(INTERVAL_VECTOR & v, INTERVAL_VECTOR & P, int const i);
RIDMME_INTERVAL_AUTODIFF LogisticFunction(INTERVAL_VECTOR & v, INTERVAL_VECTOR & P, int const i);
#endif