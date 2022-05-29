#ifndef KWADRATURY_GAUSSA_CLASS_H
#define KWADRATURY_GAUSSA_CLASS_H

class Kwadratury_Gaussa {
public:
	static double f(double x);
	static double g(double t, double alfa, double beta);
	static double kwadratura_gaussa(double* tablica_Ai, double* tablica_xi, double alfa, double beta, int ilosc_Ai);
};
#endif // !KWADRATURY_GAUSSA
