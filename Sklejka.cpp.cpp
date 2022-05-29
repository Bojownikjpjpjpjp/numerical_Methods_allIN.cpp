#include <iostream>
#include <fstream>
#include <Math.h>
#include "Gauss.h"
#include "Lagrange.h"
#include "Russel.h"
#include "LUrozwiazanie.h"
#include "Kwadratury.h"
#include "Kwadratury_Gaussa.h"

using namespace std;

int main()
{
#pragma region rozwiazanie_ukladu_rownan_gaussem
	/*
	double** A;//tworzenie tablicy 6x5 - transponowana jest
	A = new double* [6];
	for (int i = 0; i < 6; i++) {//szosta kolumna A to miejsce na wektor b
		A[i] = new double[5];
	}

	double* x;//tworzenie wektora x
	x = new double[5];
	for (int i = 0; i < 5; i++) {
		x[i] = 0;
	}
	double* b;//tworzenie wektora b
	b = new double[5];

	Gauss::fill_b(b);//sztywne wypelnienie tablicy i wektora b
	Gauss::fill_A(A, b);//tablica A ma jako ostatnią kolumnę wektor b

	//do tego miejsca dziala
	std::cout << "kom1\n";

	Gauss::wypisz(A);

	Gauss::trnasponuj_a(A);

	for (int i = 0; i < 4; i++) {
		Gauss::zeroj_kolumne(A, i);
	}

	std::cout << "kom2\n";

	Gauss::wypisz(A);
	std::cout << "kom4\n";
	Gauss::wypisz_x(x);
	Gauss::znadz_x(A, x);
	Gauss::wypisz_x(x);
	*/
#pragma endregion
#pragma region interpolacja_lagrangea
	/*
	//z zastępuje x. x jest jako tablica xów
	double z;
	int n;
	double bufor;
	double* x;


	std::cout << "podaj n: ";
	std::cin >> n;
	std::cout << "podaj x: ";
	std::cin >> z;
	std::cout << "podaj n razy xi: ";
	x = new double[5];
	for (int i = 0; i < n; i++) {
		std::cin >> bufor;
		x[i] = bufor;
	}
	for (int i = 0; i < n; i++) {
		std::cout << "\n " << x[i];
	}

	std::cout << "\n\n Ln: \n" << Lagrange::policz_Ln(z, x, n);
	return 0;
	*/
#pragma endregion
#pragma region attractor_Russela
	/*
	//parametry a,b,c
	double a = 0.2;
	double b = 0.2;
	double c = 5.7;
	//warunki początkowe 0,0,0
	double x = 0;
	double y = 0;
	double z = 0;
	double xb, yb, zb;
	//przyrost czasu
	double dt = 0.01;
	//10 000 iteracji
	int iterator = 10000;

	ofstream zapis("wyniki.txt");

	zapis << "0,0,0";
	for (int i = 0; i < iterator; i++) {
		xb = Russel::x_prim(dt, x, y, z);
		yb = Russel::y_prim(dt, x, a, y);
		zb = Russel::z_prim(dt, b, z, x, c);

		x = xb;
		y = yb;
		z = zb;
		//zapis<<"\nx: " << x << " y: " << y << " z: " << z;
		zapis << "\n" << x << "," << y << "," << z;
	}
	zapis.close();
	*/
#pragma endregion
#pragma region LU_rozwiazanie_ukladu
/*
double** L;
double** U;

//tworzenie i wypełnienie L i U
//obie macierze mają kolumnyXwiersze na początku, dlatego są potem transponowane
L = new double* [6];
for (int i = 0; i < 6; i++) {
	L[i] = new double[5];
}
U = new double* [6];
for (int i = 0; i < 6; i++) {
	U[i] = new double[5];
}

double* b;//tworzenie wektora b
b = new double[5];
LUrozwiazanie::fill_b(b);//sztywne wypelnienie tablicy i wektora b

double* x;//tworzenie wektora x i jego zerowanie
x = new double[5];
for (int i = 0; i < 5; i++) {
	x[i] = 0;
}
double* y;//tworzenie wektora y i jego zerowanie
y = new double[5];
for (int i = 0; i < 5; i++) {
	y[i] = 0;
}
LUrozwiazanie::fill_L(L, b);
LUrozwiazanie::transponuj(L);
LUrozwiazanie::wypisz(L);
LUrozwiazanie::znadz_y(L, b);
for (int i = 0; i < 5; i++) {
	y[i] = L[5][i];
}
LUrozwiazanie::wypisz_x(b);
LUrozwiazanie::wypisz(L);
LUrozwiazanie::wypisz_x(y);
LUrozwiazanie::fill_U(U, y);//recykling poprzedniego kodu z gausem
LUrozwiazanie::transponuj(U);
LUrozwiazanie::wypisz(U);
LUrozwiazanie::znadz_x(U, x);
LUrozwiazanie::wypisz(U);
LUrozwiazanie::wypisz_x(x);
*/
#pragma endregion
#pragma region Kwadratury_proste
/*
	//przedzial a, przedzial b, ilosc podzialow
	cout << Kwadratury::licz_calke(-2, 2, 566665);
	return 0;
*/
#pragma endregion
#pragma region Kwadratury_Gaussa
/*
//wagi
#pragma region tablica_Ai
double tablica_Ai[14];
tablica_Ai[0] = 1.0;
tablica_Ai[1] = 1.0;

tablica_Ai[2] = 5.0 / 9.0;
tablica_Ai[3] = 8.0 / 9.0;
tablica_Ai[4] = 5.0 / 9.0;

tablica_Ai[5] = 0.347855;
tablica_Ai[6] = 0.652145;
tablica_Ai[7] = 0.652145;
tablica_Ai[8] = 0.347855;

tablica_Ai[9] = 0.236927;
tablica_Ai[10] = 0.478629;
tablica_Ai[11] = 0.568889;
tablica_Ai[12] = 0.478629;
tablica_Ai[13] = 0.236927;
#pragma endregion
//węzły
#pragma region tablica_xi
double tablica_xi[14];
tablica_xi[0] = -0.577350;
tablica_xi[1] = 0.577350;

tablica_xi[2] = -0.774597;
tablica_xi[3] = 0.0;
tablica_xi[4] = 0.774597;

tablica_xi[5] = -0.861136;
tablica_xi[6] = -0.339981;
tablica_xi[7] = 0.339981;
tablica_xi[8] = 0.861136;

tablica_xi[9] = -0.906180;
tablica_xi[10] = -0.538469;
tablica_xi[11] = 0.0;
tablica_xi[12] = 0.538469;
tablica_xi[13] = 0.906180;
#pragma endregion

double a = -1.0;
double b = 1.0;
//f1 input
//double c = 4.5;
//double d = 0.0;
//f2 input
double c = -2.0;
double d = 2.0;

double alfa = ((b - a) / (d - c));
double beta = (((a * d) - (b * c)) / (d - c));
double wynik;
//ilosc AI ={2,3,4,5} dziala do 4 
wynik = Kwadratury_Gaussa::kwadratura_gaussa(tablica_Ai, tablica_xi, alfa, beta, 4);
cout << wynik << endl;
*/
#pragma endregion


}
