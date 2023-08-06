#include "AutomaticDifferentiation.h"

void AutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                   INTERVAL_VECTOR &,
                                                   int const),
                   Array X_inf, Array X_sup, Array P_inf, Array P_sup,
                   Matrix &GrdInf, Matrix &GrdSup, Matrix &GrdInfP,
                   Matrix &GrdSupP) {
  int n_x = X_inf.size();
  INTERVAL_VECTOR x(n_x);
  INTERVAL_VECTOR Grd(n_x);

  GrdInf.resize(n_x), GrdSup.resize(n_x);
  // Save inf and sup state intervals into INTERVAL
  for (int i = 0; i < n_x; i++) {
    x(i + 1) = INTERVAL(X_inf[i], X_sup[i]);
  }

  int n_p = P_inf.size();
  INTERVAL_VECTOR p(n_p);
  INTERVAL_VECTOR GrdP(n_x);

  if (n_p != 0) {
    GrdInfP.resize(n_x), GrdSupP.resize(n_x);
    // Save inf and sup state intervals into INTERVAL
    for (int i = 0; i < n_p; i++) {
      p(i + 1) = INTERVAL(P_inf[i], P_sup[i]);
    }
    // Calculation of inf and sup Jacobian values
    for (int i = 0; i < n_x; i++) {
      for (int j = 1; j <= n_p; j++) {
        GrdP(j) = GradientValue(fnc(x, p, i))(n_x + j);
      }
      GrdInfP[i].resize(n_p), GrdSupP[i].resize(n_p);
      for (int j = 0; j < n_p; j++) {
        GrdInfP[i][j] = Inf(GrdP(j + 1));
        GrdSupP[i][j] = Sup(GrdP(j + 1));
      }
    }
  } else {
    INTERVAL_VECTOR p(n_x);
    Matrix zeros(n_x, Array(1, 0));
    GrdInfP = zeros, GrdSupP = zeros;
  }

  for (int i = 0; i < n_x; i++) {
    for (int j = 1; j <= n_x; j++) {
      Grd(j) = GradientValue(fnc(x, p, i))(j);
    }
    GrdInf[i].resize(n_x), GrdSup[i].resize(n_x);
    for (int j = 0; j < n_x; j++) {
      GrdInf[i][j] = Inf(Grd(j + 1));
      GrdSup[i][j] = Sup(Grd(j + 1));
    }
  }
}

void MultipleAutoDiffModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                           INTERVAL_VECTOR &,
                                                           int const),
                           Matrix X_sub, MatrixInt allPerms, ArrayInt N,
                           Matrix &Jx_inf, Matrix &Jx_sup, Matrix &Jp_inf,
                           Matrix &Jp_sup, ArrayBool &AnyNanElement) {

  int NumPers = allPerms.size();
  int n_x = Jx_inf[0].size();
  int n_p = Jp_inf[0].size();

  MatrixList GrdInfList(NumPers), GrdSupList(NumPers);
  MatrixList GrdInfPList(NumPers), GrdSupPList(NumPers);

  N = getCumulativeSum(N);

  bool ContinueFlag = true;

#pragma omp parallel for
  for (int i = 0; i < NumPers; i++) {

    Array X_inf(n_x);
    Array X_sup(n_x);

    Array P_inf(n_p);
    Array P_sup(n_p);

    Matrix GrdInf(n_x, Array(n_x, 0));
    Matrix GrdSup(n_x, Array(n_x, 0));

    Matrix GrdInfP(n_x, Array(n_p, 0));
    Matrix GrdSupP(n_x, Array(n_p, 0));

    if (ContinueFlag == true) {
      ArrayInt row = allPerms[i];
      for (int j = 0; j < n_x + n_p; j++) {
        int idx = row[j];
        if (j < n_x) {
          X_inf[j] = X_sub[N[j] + idx - 1][0];
          X_sup[j] = X_sub[N[j] + idx - 1][1];
        } else {
          P_inf[j - n_x] = X_sub[N[j] + idx - 1][0];
          P_sup[j - n_x] = X_sub[N[j] + idx - 1][1];
        }
      }
      AutoDiffModel(fnc, X_inf, X_sup, P_inf, P_sup, GrdInf, GrdSup, GrdInfP,
                    GrdSupP);

      // check any element of gradient/jacobian are Nan
      ArrayBool GrdNan = checkForNaN(GrdInf, GrdSup, 1);
      ArrayBool GrdPNan = checkForNaN(GrdInfP, GrdSupP, 1);

      if (any_of(GrdNan.begin(), GrdNan.end(),
                 [](bool value) { return value; })) {
        cout << "Identified: jacobian or gradient element in NAN" << endl;
        AnyNanElement = GrdNan;
        ContinueFlag = false;
      } else if (any_of(GrdPNan.begin(), GrdPNan.end(),
                        [](bool value) { return value; })) {
        cout << "Identified: jacobian or gradient element in NAN" << endl;
        AnyNanElement = GrdPNan;
        ContinueFlag = false;
      }
      GrdInfList[i] = GrdInf;
      GrdSupList[i] = GrdSup;
      GrdInfPList[i] = GrdInfP;
      GrdSupPList[i] = GrdSupP;
    }
  }
  if (ContinueFlag == false) {
    Jx_inf = Matrix(n_x, Array(n_x, -numeric_limits<double>::infinity()));
    Jx_sup = Matrix(n_x, Array(n_x, numeric_limits<double>::infinity()));
    Jp_inf = Matrix(n_x, Array(n_p, -numeric_limits<double>::infinity()));
    Jp_sup = Matrix(n_x, Array(n_p, numeric_limits<double>::infinity()));
    return;
  }
  Jx_inf = findMinMaxMatrix(GrdInfList, 0);
  Jx_sup = findMinMaxMatrix(GrdSupList, 1);

  Jp_inf = findMinMaxMatrix(GrdInfPList, 0);
  Jp_sup = findMinMaxMatrix(GrdSupPList, 1);

  cout << "" << endl;
}

void ToDoubleEvaluationModel(RIDMME_INTERVAL_AUTODIFF (*fnc)(INTERVAL_VECTOR &,
                                                             INTERVAL_VECTOR &,
                                                             int const),
                             double t, Array &x, Array &p) {
  int n_x = x.size();
  INTERVAL_VECTOR X(n_x);
  for (int i = 0; i < n_x; i++) {
    X(i + 1) = INTERVAL(x[i], x[i]);
  }

  int n_p = p.size();
  INTERVAL_VECTOR P(n_p);
  for (int i = 0; i < n_p; i++) {
    P(i + 1) = INTERVAL(p[i], p[i]);
  }

  Array Sol(n_x);
  for (int i = 0; i < n_x; i++) {
    Sol[i] = Inf(FunctionValue(fnc(X, P, i)));
  }
  x = Sol;
}

RIDMME_INTERVAL_AUTODIFF MergeTwoIntervalArr(INTERVAL_VECTOR &v,
                                             INTERVAL_VECTOR &P) {
  auto FactorsSize = Dimension(v) + Dimension(P);
  INTERVAL_VECTOR Sol(FactorsSize);

  for (int i = 0; i < Dimension(v); i++) {
    Sol(i + 1) = v(i + 1);
  }
  for (int i = 0; i < Dimension(P); i++) {
    Sol(i + Dimension(v) + 1) = P(i + 1);
  }
  return Sol;
}