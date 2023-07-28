#include "Metrics.h"

double MetricEval(array ydata, array Constraint, int ErrorMetric){
    double ydata_low = ydata[0];
    double ydata_up = ydata[1];
    double Constraint_low = Constraint[0];
    double Constraint_up = Constraint[1];
    double Error;

    switch (ErrorMetric) {
        case 1: // Hausdorff
            Error = max(abs(Constraint_low - ydata_low), abs(Constraint_up - ydata_up));
            break;
        case 2: // CoreBased
            Error = 1 / sqrt(2) * sqrt(pow(ydata_low - Constraint_low, 2) + pow(ydata_up - Constraint_up, 2));
            break;
        case 3: // Bertoluzza
            double C_c = (Constraint_low + Constraint_up) / 2;
            double X_c = (ydata_low + ydata_up) / 2;
            double C_r = C_c - Constraint_low;
            double X_r = X_c - ydata_low;
            Error = sqrt(pow(X_c - C_c, 2) + (1.0 / 3.0) * pow(X_r - C_r, 2));
    }
    return Error;
}

array TimeSeriesError(array (*fnc)(array X_inf, array X_sup), matrix X_inf, matrix X_sup, matrix Constraint, int ErrorMetric){
    int t_f = X_inf.size();
    array ydata;
    array Error(t_f,0);
    for(int i = 0; i < t_f; i++){
        ydata = fnc(X_inf[i], X_sup[i]);
        Error[i] = MetricEval(ydata, Constraint[i], ErrorMetric);
    }
    return Error;
}

void MultipleError(array (*fnc)(array X_inf, array X_sup), matrixList X_inf, matrixList X_sup, matrix Constraint, int ErrorMetric, string Name){
    int n_l = X_inf.size();
    array Error;

    for(int i = 0; i < n_l; i++){
        Error = TimeSeriesError(fnc, X_inf[i], X_sup[i], Constraint, ErrorMetric);
        saveArrayToCSV(Error, array(), "Error/Error" + Name + to_string(i) + ".csv");
    }
}