#include <math.h>
#include "Kwadratury.h"
double Kwadratury::f(double x) {
	//return ((x*x)*pow(sin(x),3));
	return (exp((x * x)) * (x - 1));

}
double Kwadratury::licz_calke(double a, double b, int n) {
	double suma = 0;
	double a1 = a;
	double dx = ((b - a) / n);
	if (dx < 0) {
		dx *= -1;
	}
	double b1 = a1 + dx;
	double suma_podstaw = 0;
	for (int i = 0; i < n; i++) {
		suma_podstaw = f(a1) + f(b1);
		suma += (suma_podstaw * dx / 2);
		a1 += dx;
		b1 += dx;
	}
	return suma;
}