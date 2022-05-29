#include "LUrozwiazanie.h"

void LUrozwiazanie::fill_b(double* b) {
	b[0] = 74.64;
	b[1] = -40.26;
	b[2] = -2.32;
	b[3] = 12.6;
	b[4] = -8.9;
}
void LUrozwiazanie::fill_A(double** A, double* b) {
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
void LUrozwiazanie::fill_L(double** A, double* b) {//A pozosta³o bo mi sie nie chcia³o zmieniaæ
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			A[j][i] = 0;
		}
	}
	A[0][0] = 1;

	A[1][0] = -3;
	A[1][1] = 1;

	A[2][0] = 2;
	A[2][1] = -0.7142857;
	A[2][2] = 1;


	A[3][0] = -6;
	A[3][1] = 2.2857143;
	A[3][2] = -1.042918;
	A[3][3] = 1;

	A[4][0] = 5;
	A[4][1] = -1.8571429;
	A[4][2] = 1.551502;
	A[4][3] = -1.350949;
	A[4][4] = 1;

	A[5][0] = b[0];
	A[5][1] = b[1];
	A[5][2] = b[2];
	A[5][3] = b[3];
	A[5][4] = b[4];

}
void LUrozwiazanie::fill_U(double** A, double* y) {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			A[j][i] = 0;
		}
	}
	A[0][0] = 1;
	A[0][1] = -3;
	A[0][2] = 4;
	A[0][3] = 6.8;
	A[0][4] = 9;

	A[1][1] = -7;
	A[1][2] = 16.6;
	A[1][3] = 26.7;
	A[1][4] = 17;

	A[2][2] = 6.657143;
	A[2][3] = -2.928571;
	A[2][4] = -10.857143;

	A[3][3] = -23.782833;
	A[3][4] = 2.919742;

	A[4][4] = 3.360733;

	A[5][0] = y[0];
	A[5][1] = y[1];
	A[5][2] = y[2];
	A[5][3] = y[3];
	A[5][4] = y[4];
}
void LUrozwiazanie::wypisz(double** A) {
	std::cout << "\n";
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 6; j++) {
			std::cout << A[j][i] << " ";
		}
		std::cout << "\n";
	}
}
void LUrozwiazanie::wypisz_x(double* x) {
	std::cout << "\n";
	for (int i = 0; i < 5; i++) {
		std::cout << x[i] << " ";
	}
	std::cout << "\n";
}
void LUrozwiazanie::zeroj_kolumne(double** A, int iteracja) {//najpier kol = 0
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
void LUrozwiazanie::znadz_x(double** A, double* x) {
	x[4] = A[5][4] / A[4][4];//[kolumna][wiersz], A[6][y] =b
	x[3] = (A[5][3] - (x[4] * A[4][3])) / A[3][3];
	x[2] = (A[5][2] - (x[3] * A[3][2]) - (x[4] * A[4][2])) / A[2][2];
	x[1] = (A[5][1] - (x[2] * A[2][1]) - (x[3] * A[3][1]) - (x[4] * A[4][1])) / A[1][1];
	x[0] = (A[5][0] - (x[1] * A[1][0]) - (x[2] * A[2][0]) - (x[3] * A[3][0]) - (x[4] * A[4][0])) / A[0][0];
}
void LUrozwiazanie::znadz_y(double** A, double* b) {
	double suma;//dla macierzy L
	A[5][0] = b[0];
	for (int i = 1; i < 5; i++) {//macieŸ jest dalej [kolumna][wiersz]
		suma = 0;
		for (int j = 0; j < i; j++) {
			//cout << "A j i  " << A[j][i] << "    " << "A 5 j  " << A[5][j];
			suma = suma + (A[j][i] * A[5][j]);
		}
		//cout<<"suma: " << suma << "\n";
		A[5][i] = b[i] - suma;
	}
}
void LUrozwiazanie::transponuj(double** A) {
	double transpose[5][5];
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j) {
			transpose[j][i] = A[i][j];
		}
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j) {
			A[j][i] = transpose[j][i];
		}

}