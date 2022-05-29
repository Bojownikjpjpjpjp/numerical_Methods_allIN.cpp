#ifndef RUSSEL_CLASS_H
#define RUSSEL_CLASS_H

#include <iostream>
#include <Math.h>
class Russel {
public:
    static double x_prim(double dt, double x, double y, double z);
    static double y_prim(double dt, double x, double a, double y);
    static double z_prim(double dt, double b, double z, double x, double c);

};
#endif // !RUSSEL_CLASS_H

