#include <iostream>
#include <Math.h>
#include "Gauss.h"



void Gauss::fill_b(double* b) {
	b[0] = 74.64;
	b[1] = -40.26;
	b[2] = -2.32;
	b[3] = 12.6;
	b[4] = -8.9;
}
void Gauss::fill_A(double** A, double* b) {
	A[0][0] = 1;
	A[0][1] = -3;
	A[0][2] = 4;
	A[0][3] = 6.8;
	A[0][4] = 9;

	A[1][0] = -3;
	A[1][1] = 2;
	A[1][2] = 4.6;
	A[1][3] = 6.3;
	A[1][4] = -10;

	A[2][0] = 2;
	A[2][1] = -1;
	A[2][2] = 2.8;
	A[2][3] = -8.4;
	A[2][4] = -5;

	A[3][0] = -6;
	A[3][1] = 2;
	A[3][2] = 7;
	A[3][3] = -0.5;
	A[3][4] = -0.9;

	A[4][0] = 5;
	A[4][1] = -2;
	A[4][2] = -0.5;
	A[4][3] = 12;
	A[4][4] = -4;

	A[5][0] = b[0];
	A[5][1] = b[1];
	A[5][2] = b[2];
	A[5][3] = b[3];
	A[5][4] = b[4];

}

void Gauss::wypisz(double** A) {
	std::cout << "\n";
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 6; j++) {
			std::cout << A[j][i] << " ";
		}
		std::cout << "\n";
	}
}
void Gauss::wypisz_x(double* x) {
	std::cout << "\n";
	for (int i = 0; i < 5; i++) {
		std::cout << x[i] << " ";
	}
	std::cout << "\n";
}
void Gauss::zeroj_kolumne(double** A, int iteracja) {//najpier kol = 0
	double iloraz;
	//	double czynnik2;
	wypisz(A);
	std::cout << "kom3\n";

	if (iteracja <= 3) {
		for (int j = iteracja; j < 4; j++) {//wiersz
			iloraz = A[iteracja][j + 1] / A[iteracja][iteracja];
			for (int i = 0; i < 6; i++) {//kolumna
				A[i][j + 1] -= iloraz * A[i][iteracja];
			}
		}
	}
}
void Gauss::znadz_x(double** A, double* x) {
	x[4] = A[5][4] / A[4][4];//[kolumna][wiersz], A[6][y] =b
	x[3] = (A[5][3] - (x[4] * A[4][3])) / A[3][3];
	x[2] = (A[5][2] - (x[3] * A[3][2]) - (x[4] * A[4][2])) / A[2][2];
	x[1] = (A[5][1] - (x[2] * A[2][1]) - (x[3] * A[3][1]) - (x[4] * A[4][1])) / A[1][1];
	x[0] = (A[5][0] - (x[1] * A[1][0]) - (x[2] * A[2][0]) - (x[3] * A[3][0]) - (x[4] * A[4][0])) / A[0][0];
}
void Gauss::trnasponuj_a(double** A) {
	double transpose[5][5];
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j) {
			transpose[j][i] = A[i][j];
		}
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j) {
			A[j][i] = transpose[j][i];
		}

	//	for (int wiersze = 0; wiersze < 5; wiersze++) {
	//		for (int kolumny = 0; kolumny < 5; kolumny++) {
	//
	//		}
	//	}
}


