#include <iostream>
#include <Math.h>
#include "Russel.h"

double Russel::x_prim(double dt, double x, double y, double z) {
    return (x + (dt * (-y - z)));
}
double Russel::y_prim(double dt, double x, double a, double y) {
    return (y + (dt * (x + (a * y))));
}
double Russel::z_prim(double dt, double b, double z, double x, double c) {
    return z + (dt * (b + (z * (x - c))));
}



