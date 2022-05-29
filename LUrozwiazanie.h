#ifndef LUROZWIAZANIE_CLASS_H
#define LUROZWIAZANIE_CLASS_H

#include <iostream>
#include <Math.h>
class LUrozwiazanie {
public:
	static void fill_b(double* b);
	static void fill_A(double** A, double* b);
	static void fill_L(double** A, double* b);
	static void fill_U(double** A, double* y);
	static void wypisz(double** A);
	static void wypisz_x(double* x);
	static void zeroj_kolumne(double** A, int iteracja);
	static void znadz_x(double** A, double* x);
	static void znadz_y(double** A, double* b);
	static void transponuj(double** A);
};
#endif // !LUROZWIAZANIE

