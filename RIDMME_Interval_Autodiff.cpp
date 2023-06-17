#include "RIDMME_Interval_Autodiff.h"
#include "Functions.h"

#ifdef AUTODIFF_COMPUTE_HESSIAN
extern INTERVAL_MATRIX TransposedMul (CONST INTERVAL_VECTOR & a,
                                    CONST INTERVAL_VECTOR & b);
#endif

RIDMME_INTERVAL_AUTODIFF Power(CONST RIDMME_INTERVAL_AUTODIFF &g, INT p)
{
  INT dim = Dimension(g);
  RIDMME_INTERVAL_AUTODIFF h(dim, FALSE);

  h.fkt() = Power(g.fkt(), p);
  h.grd() = p * Power(g.fkt(), p-1) * g.grd();
#ifdef AUTODIFF_COMPUTE_HESSIAN
  if (g.ComputeHessian) {
    h.hessian() = p * (Power(g.fkt(), p-1) * g.hessian() + TransposedMul (g.grd(), g.grd()));
    h.ComputeHessian = TRUE;
  }
#endif
  return h;
}

RIDMME_INTERVAL_AUTODIFF Power(CONST INTERVAL_AUTODIFF& g, INT p)
{
  RIDMME_INTERVAL_AUTODIFF h = (RIDMME_INTERVAL_AUTODIFF)g;
  return Power(h, p);
}