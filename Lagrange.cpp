#include <iostream>
#include <Math.h>
#include "Lagrange.h"

double Lagrange::f(double x) {
    return  (x - 5) * x + 14;
}

double Lagrange::policz_li(double* x, double z, int i, int n) {
    double iloraz = 1;
    for (int j = 0; j < n; j++) {
        if (j == i) {
            continue;
        }
        iloraz = iloraz * ((z - x[j]) / (x[i] - x[j]));
    }
    return iloraz;
}

double Lagrange::policz_Ln(double z, double* x, int n) {
    double suma = 0;
    for (int i = 0; i < n; i++) {
        suma = suma + f(x[i]) * policz_li(x, z, i, n);
    }
    return suma;
}



