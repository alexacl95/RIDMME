#ifndef MO_H
#define MO_H
#include "Matrix.h"
#include "AutomaticDifferentiation.h"

INTERVAL_AUTODIFF ModelFunction1(INTERVAL_VECTOR & v, int const i);
INTERVAL_AUTODIFF LogisticFunction(INTERVAL_VECTOR & v, int const i);


#endif