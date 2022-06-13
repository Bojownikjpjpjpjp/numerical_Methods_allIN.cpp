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
double f(double T, double alfa, double beta) {
	return (alfa*(pow(T,4)-beta));
}

#pragma region Ortogonalizacja_funkcje
/*

double(*f[6])(double);
double(*g[6])(double, double**, int, int);
*/
#pragma region Funkcje_przestrzeni_F
double f0(double x) {
	return pow(x, 0);
}
double f1(double x) {
	return pow(x, 1);
}
double f2(double x) {
	return pow(x, 2);
}
double f3(double x) {
	return pow(x, 3);
}
double f4(double x) {
	return pow(x, 4);
}

double f5(double x) {
	return pow(x, 5);
}

#pragma endregion
#pragma region Funkcje_przestrzeni_G
double g0(double x, double** bazaG, int i_UP, int j_UP) {
	double suma = 0;
	//i_UP--;
	for (int j = 0; j < i_UP; j++) {
		suma += bazaG[i_UP][j] * pow(x, j);
	}
	return suma;
	//return pow((x * bazaG[i_UP][i_UP]), i_UP);
	//sume po j G[i][j] * pow(x, j)
}
double g1(double x, double** bazaG, int i_UP, int j_UP) {
	double suma = 0;
	i_UP--;
	//cout << x << " " << i_UP << " " << j_UP;
	for (int j = 0; j <= i_UP; j++) {
		suma += bazaG[i_UP][j] * pow(x, j);
		//cout <<j<<bazaG[i_UP][j]<<" " << suma;
	}
	
	return suma;
	//return pow((x * bazaG[i_UP][i_UP]), i_UP);
}
double g2(double x, double** bazaG, int i_UP, int j_UP) {
	double suma = 0;
	//i_UP--;
	
	for (int j = 0; j < i_UP; j++) {
		suma += bazaG[i_UP][j] * pow(x, j);
	}
	return suma;
	//return pow((x * bazaG[i_UP][i_UP]), i_UP);
}
double g3(double x, double** bazaG, int i_UP, int j_UP) {
	double suma = 0;
	//i_UP--;
	for (int j = 0; j < i_UP; j++) {
		suma += bazaG[i_UP][j] * pow(x, j);
	}
	return suma;
	//return pow((x * bazaG[i_UP][i_UP]), i_UP);
}
double g4(double x, double** bazaG, int i_UP, int j_UP) {
	double suma = 0;
	//i_UP--;
	for (int j = 0; j < i_UP; j++) {
		suma += bazaG[i_UP][j] * pow(x, j);
	}
	return suma;
	//return pow((x * bazaG[i_UP][i_UP]), i_UP);
}
double g5(double x, double** bazaG, int i_UP, int j_UP) {
	double suma = 0;
	//i_UP--;
	for (int j = 0; j < i_UP; j++) {
		suma += bazaG[i_UP][j] * pow(x, j);
	}
	return suma;
	//return pow((x * bazaG[i_UP][i_UP]), i_UP);
}

#pragma endregion

double licz_calke_z_gg(double a, double b, int n, int i_UP, int j_UP, double (*G1[])(double, double**, int, int), double (*G2[])(double, double**, int, int), double** bazaG) {
	//a,b przedział, n ilość podpodziałów
	double suma = 0;
	double a1 = a;
	double dx = ((b - a) / n);
	if (dx < 0) {
		dx *= -1;
	}
	double b1 = a1 + dx;

	double suma_podstaw = 0;
	for (int i = 0; i < n; i++) {
		//suma_podstaw = ((F1(a1)*F2(a1)) + (F1(b1)*F2(b1)));
		suma_podstaw = ((G1[1](a1, bazaG, i_UP, j_UP) * G1[1](a1, bazaG, i_UP, j_UP)) + (G1[1](b1, bazaG, i_UP, j_UP) * (G2[1](b1, bazaG, i_UP, j_UP))));

		suma += ((suma_podstaw * dx) / 2);
		a1 += dx;
		b1 += dx;
	}
	cout << "\ncalka: " << suma;
	return suma;
}

double licz_calke_z_gf(double a, double b, int n, int i_UP, int j_UP, double (*G[])(double, double**, int,int), double (*F2)(double), double** bazaG) {
	//a,b przedział, n ilość podpodziałów
	double suma = 0;
	double a1 = a;
	double dx = ((b - a) / n);
	if (dx < 0) {
		dx *= -1;
	}
	double b1 = a1 + dx;

	double suma_podstaw = 0;
	for (int i = 0; i < n; i++) {
		suma_podstaw = ((G[1](a1, bazaG, i_UP, j_UP) * F2(a1)) + (G[1](a1, bazaG, i_UP, j_UP) * F2(b1)));
		suma += ((suma_podstaw * dx) / 2);
		a1 += dx;
		b1 += dx;
	}
	cout << "\ncalka: " << suma;
	return suma;
}

double licz_wspolczynniki(double** bazaF, double** bazaG, double a, double b, double (*f[])(double), double (*g[])(double, double**, int, int), int n) {
	double licznik = 0.0;
	double mianownik = 0.0;
	double w_sumy = 0.0;
	double suma = 0.0;

	for (int i = 1; i < 6; i++) {//pętla liczenia współczynników
		
		for (int j = 0; j < i; j++) {//pętla sumy
			licznik = licz_calke_z_gf(a, b, n, i, j, g, f[i], bazaG);
			cout << "\n Licznik: " << licznik;
			mianownik = licz_calke_z_gg(a, b, n, i, j, g, g, bazaG);
			cout << "\n mianownik: " << mianownik;
			w_sumy = ((licznik / mianownik)* bazaG[j][j]);
			suma += w_sumy;
			cout << "\n wyraz sumy: " << w_sumy;
			cout << "\n suma: " << suma;
			//cout << "\n wyraz bazy G: " << bazaG[i][i];
			//cout << "\n wyraz bazy F: " << bazaF[i][i];
			bazaG[i][j] = w_sumy;
			w_sumy = 0.0;
			//bazaG[i][j] = w_sumy;
		}

		bazaG[i][i] = (bazaF[i][i] - suma);
		cout << "\n po przypisaniu wyraz bazy G: " << bazaG[i][i];
		cout << "\n ";
		suma = 0.0;
	}
	return suma;
}

#pragma endregion
int main()
{
#pragma region rownania_rozniczkowe
	// dT/dt=alfa(T^4-beta) //y'=f(x,y)
	// T(t)= alfa(T^4-beta)=f(y)
	//T(0)=T0=1200K
	// T po 300 sekundach
	double T = 1200;
	double t = 0;
	double h = 0.1;
	double alfa = -1e-12;
	double beta = 0;
	ofstream zapis("wyniki_rr.txt");
	zapis << "czas,temperatura";
	zapis << "0,0";
	for (t; t < 300; t += h) {
		T = (T+(h*(f(T, alfa, beta))));
		zapis << "\n" << t << "," << T;
		cout << "\nT = " << T;
	}

	cout << "\nTk = " << T;

	zapis.close();
#pragma endregion


#pragma region Ortogonalizacje
/*
#pragma region tworzenie_i_zerowanie_baz
	
	//baza funkcji. całość w postaci wielomianu :a0*1, a1*x, a2*x^2...a4*x^4
	//tworzenie bazy
	double** bazaF;
	double** bazaG;
	double** bazaG2;
	bazaF = new double* [6];
	for (int i = 0; i < 6; i++) {
		bazaF[i] = new double[6];
	}
	//zerowanie bazy
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			bazaF[i][j] = 0.0;
		}
	}
	bazaG = new double* [6];
	for (int i = 0; i < 6; i++) {
		bazaG[i] = new double[6];
	}
	//zerowanie bazy
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			bazaG[i][j] = 0.0;
		}
	}
	bazaG2 = new double* [6];
	for (int i = 0; i < 6; i++) {
		bazaG2[i] = new double[6];
	}
	//zerowanie bazy
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			bazaG2[i][j] = 0.0;
		}
	}
#pragma endregion
#pragma region przypisanie_wskaźników_funkcji
	//funkcje przestrzeni: 1,x,x^2,x^3,x^4
	f[0] = f0;
	f[1] = f1;
	f[2] = f2;
	f[3] = f3;
	f[4] = f4;
	f[5] = f5;
	//funkcje przestrzeni: 1,x,x^2,x^3,x^4
	g[0] = g0;
	g[1] = g1;
	g[2] = g2;
	g[3] = g3;
	g[4] = g4;
	g[5] = g5;
#pragma endregion
#pragma region Ustawienie_Diagonalnej_bazy_standardowej_i_pierwszego_diagonalnego_wyrazu_macierzy_g
//baza F to baza podstawowa 1,x,x^2,x^3,x^4
//baza g ma w sobie współczynniki ortogonalizacji bazy f, i po algorytmie jest bazą ortogonalną przestrzeni funkcyjnej

//[wektory macierzy][współczynniki]
//współczynniki bazy standardowej (na diagonalnej 1, reszta 0
	bazaF[0][0] = 1;
	bazaF[1][1] = 1;
	bazaF[2][2] = 1;
	bazaF[3][3] = 1;
	bazaF[4][4] = 1;
	bazaF[5][5] = 1;
	//[wektory macierzy][współczynniki]
	//współczynniki
	bazaG[0][0] = bazaF[0][0];

#pragma endregion
	double a = -1.0;
	double b = 1.0;
	
	cout<<licz_wspolczynniki(bazaF, bazaG, a, b, f, g , 20000);
	

#pragma region wypisanie_baz
	cout << "baza Podstawowa: \n";
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			//if (bazaF[i][j] != 0) {
			//	cout << bazaF[i][j] << "\n";
				cout << bazaF[i][j] << " ";
			//}
			//else {
			//	continue;
			//}
		}
		cout << " \n";
	}
	//cout << "\nG [1][1] "<<bazaG[1][1]<<"\n";
	cout << "\nbaza ortogonalna: \n";
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			//if (bazaG[i][j] != 0) {
				//cout << bazaG[i][j] << "\n";
				cout << bazaG[i][j] << " ";
			//}
			//else {
			//	continue;
			//}
		}
		cout <<" \n";
	}
#pragma endregion
*/
#pragma region Reortogonalizacja_próba
/*
	bazaG2[0][0] = bazaG[0][0];
	bazaG2[1][1] = bazaG[1][1];
	cout << licz_wspolczynniki(bazaG, bazaG2, a, b, f, 9000);


#pragma region wypisanie_baz_2
	cout << "baza Podstawowa: \n";
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			//if (bazaF[i][j] != 0) {
			//	cout << bazaF[i][j] << "\n";
			cout << bazaG[i][j] << " ";
			//}
			//else {
			//	continue;
			//}
		}
		cout << " \n";
	}
	//cout << "\nG [1][1] "<<bazaG[1][1]<<"\n";
	cout << "\nbaza ortogonalna: \n";
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			//if (bazaG[i][j] != 0) {
				//cout << bazaG[i][j] << "\n";
			cout << bazaG2[i][j] << " ";
			//}
			//else {
			//	continue;
			//}
		}
		cout << " \n";
	}
#pragma endregion
*/
#pragma endregion

#pragma endregion
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

return 0;
}
