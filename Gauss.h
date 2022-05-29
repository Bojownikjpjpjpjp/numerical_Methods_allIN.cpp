#ifndef GAUSS_CLASS_H
#define GAUSS_CLASS_H

#include <iostream>
#include <Math.h>
class Gauss {
public:
	static void fill_b(double* b);
	static void fill_A(double** A, double* b);
	static void wypisz(double** A);
	static void wypisz_x(double* x);
	static void zeroj_kolumne(double** A, int iteracja);
	static void znadz_x(double** A, double* x);
	static void trnasponuj_a(double** A);
};
#endif // !GAUSS_H

