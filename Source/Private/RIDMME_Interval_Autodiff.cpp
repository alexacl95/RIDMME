#include "RIDMME_Interval_Autodiff.h"
#include "Functions.h"

#ifdef AUTODIFF_COMPUTE_HESSIAN
extern INTERVAL_MATRIX TransposedMul (CONST INTERVAL_VECTOR & a,
                                    CONST INTERVAL_VECTOR & b);
#endif

//Log function
RIDMME_INTERVAL_AUTODIFF Log(CONST RIDMME_INTERVAL_AUTODIFF &g){
  INT dim = Dimension(g);
  RIDMME_INTERVAL_AUTODIFF h(dim, FALSE);

  h.fkt() = Log(g.fkt());
  h.grd() = g.grd()/g.fkt();
#ifdef AUTODIFF_COMPUTE_HESSIAN
  if (g.ComputeHessian) {
    h.hessian() = TransposedMul(g.grd(), g.grd());
    h.ComputeHessian = TRUE;
  }
#endif
  return h;
}

RIDMME_INTERVAL_AUTODIFF Log(CONST INTERVAL_AUTODIFF& g){
  RIDMME_INTERVAL_AUTODIFF h = (RIDMME_INTERVAL_AUTODIFF)g;
  return Log(h);
}


//Power function of real n
RIDMME_INTERVAL_AUTODIFF Power(CONST RIDMME_INTERVAL_AUTODIFF &g, INT p){
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

RIDMME_INTERVAL_AUTODIFF Power(CONST INTERVAL_AUTODIFF& g, INT p){
  RIDMME_INTERVAL_AUTODIFF h = (RIDMME_INTERVAL_AUTODIFF)g;
  return Power(h, p);
}

//Power of f(x)^g(x)
RIDMME_INTERVAL_AUTODIFF Power(CONST RIDMME_INTERVAL_AUTODIFF &g, CONST RIDMME_INTERVAL_AUTODIFF &f){
  INT dim = Dimension(g);
  RIDMME_INTERVAL_AUTODIFF h(dim, FALSE);
  
  h.fkt() = Power(g.fkt(), f.fkt());
  h.grd() = Power(g.fkt(), f.fkt()-1)*(Log(g.fkt())*g.fkt()*f.grd() + f.fkt()*g.grd());
#ifdef AUTODIFF_COMPUTE_HESSIAN
  if (g.ComputeHessian) {
    h.hessian() = Power(g.fkt(), 1) * g.hessian() + TransposedMul (g.grd(), g.grd());
    h.ComputeHessian = TRUE;
  }
#endif
  return h;
}

RIDMME_INTERVAL_AUTODIFF Power(CONST INTERVAL_AUTODIFF& g, CONST RIDMME_INTERVAL_AUTODIFF &f){
  RIDMME_INTERVAL_AUTODIFF h = (RIDMME_INTERVAL_AUTODIFF)g;
  RIDMME_INTERVAL_AUTODIFF h1 = (RIDMME_INTERVAL_AUTODIFF)f;
  bool conditionBase = Inf(h.fkt()) < 0 ||  Sup(h.fkt()) < 0;
  if(conditionBase){
    RIDMME_INTERVAL_AUTODIFF IntervalNan = h;
    IntervalNan.fkt() = INTERVAL(NAN, NAN);
    IntervalNan.ComputeHessian = TRUE;
    cout << "Negative base in power function" << endl;
    return IntervalNan;
  }
  return Power(h, h1);
}
