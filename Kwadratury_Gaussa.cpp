#include "Kwadratury_Gaussa.h"
#include <math.h>
#include <iostream>
double Kwadratury_Gaussa::f(double x) {//funkcja 1
	//return ((x * x)*(pow(sin(x), 3)));//f1
	return (exp((x * x)) * (x - 1));//f2
}
double Kwadratury_Gaussa::g(double t, double alfa, double beta) {//funkcja 3
	//double wynik = (f((alfa * t) + beta));
	double wynik = (f(t)); //dzia³a dla exp		f2
	//double wynik = f(alfa * t);//dzia³a lepiej dla x^2 *sin^3		f1
	return wynik;
}

double Kwadratury_Gaussa::kwadratura_gaussa(double* tablica_Ai, double* tablica_xi, double alfa, double beta, int ilosc_Ai) {
	//tablica A - wagi
	//tablica x - wêz³y
	double suma = 0.0;
	double calka;
	int setter;//wskazuje miejsce w tablicy AI
	int iteracje;//okresla ile razy ma byæ iterowana funkcja (ile sk³adników)
#pragma region ustawienie_poczatku_tablicy_i_ilosci_obrotow_petli
	if (ilosc_Ai == 2) {
		setter = 0;
		iteracje = 2;
	}
	else if (ilosc_Ai == 3) {
		setter = 2;
		iteracje = 3;
	}
	else if (ilosc_Ai == 4) {
		setter = 5;
		iteracje = 4;
	}
	else if (ilosc_Ai == 5) {
		setter = 9;
		iteracje = 5;
	}
	else {
		std::cout << "!!!cos sie zepsulo!!!";
		return 0.0;
	}
#pragma endregion

	double Wagi[13];//Bi = tablica_Ai[i] / alfa;
	double Wezly[13];
	//przypisanie i policzenie Ti
	for (int i = setter; i < (setter + iteracje); i++) {
		Wagi[i] = tablica_Ai[i] / alfa;
		//Wagi[i] = tablica_Ai[i];
	}
	for (int i = setter; i < (setter + iteracje); i++) {
		Wezly[i] = ((tablica_xi[i] - beta) / alfa);
		//Wezly[i] = (tablica_xi[i] / alfa);
	}

	for (int i = setter; i < (setter + iteracje); i++) {
		calka = (Wagi[i] * g(Wezly[i], alfa, beta));
		//calka = Bi[i] * (v(Ti[i], alfa, beta)+ g(Ti[i], alfa, beta));

		//cout << "\n iteracja nr: " << i;

		suma += calka;
		std::cout << "\n suma+calka: " << suma << "\n";
	}

	//std::cout << "\ntu doszlo\n";
	return suma;
}