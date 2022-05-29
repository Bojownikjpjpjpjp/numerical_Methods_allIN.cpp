#ifndef LAGRANGE_CLASS_H
#define LAGRANGE_CLASS_H

class Lagrange {
public:
    static double f(double x);
    static double policz_li(double* x, double z, int i, int n);
    static double policz_Ln(double z, double* x, int n);
};
#endif // !LAGRANGE_CLASS_H
