#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <string>
#include <math.h>
#include <fstream>

using namespace std;

double TrigonometricInterpol(double t, double* X, double* g, int n)
{
	//числовые коэффициенты
	double A0 = 0;
	double Ak = 0;
	double Bk = 0;

	double sum = 0;
	double result = 0;

	for (int i = 0; i < 2 * n + 1; i++){

		sum += g[i];
	}

	A0 = sum / (2 * n + 1);
	result += A0;

	for (int k = 1; k <= n; k++){
			sum = 0;

			for (int i = 0; i < 2 * n + 1; i++) {
				sum += g[i] * cos(k * X[i]);
			}
			Ak = 2 * sum / (2 * n + 1);
			sum = 0;

			for (int i = 0; i < 2 * n + 1; i++) {
				sum += g[i] * sin(k * X[i]);
			}

			Bk = 2 * sum / (2 * n + 1);
			result += Ak * cos(k * t) + Bk * sin(k * t);

	}
	
		return result;
}


double G_t(double t)
{
	return atan(t)/(1+t*t);//функция
}

int main() {

	double a = 0;
	double b = 2;
	int c = 1000;

	ofstream fout;
	ofstream fout2;

	fout.open("delta.txt");//погрешность дельта
	fout2.open("graph.txt");//расчет интерполяции

	for (int n = 1; n < 300; n++) {

		double h = (b - a) / n;//сетка

		double* t = new double[2 * n + 1];
		double* g = new double[2 * n + 1];

		for (int i = 0; i < 2 * n + 1; i++) {
			t[i] = 2 * M_PI * (i) / (2 * n + 1);//узлы

			g[i] = G_t(t[i]);//значения узлов
		}

		double* t_rep = new double[c];//замена

		double step = (b - a) / c;

		for (int i = 0; i < c - 1; i++){

			t_rep[i] = a + step * i;//рассматриваемые участки интерполяции
		}
			double MaxDelta = 0;

			for (int i = 0; i < c - 1; i++) {

				double delta = abs(TrigonometricInterpol(t_rep[i], t, g, n) - G_t(t_rep[i]));//дельта

				if (n == 100) {

					fout2 << delta << endl;
				}

				if (delta > MaxDelta) {

					MaxDelta = delta;
				}

			}

			cout << n % 10;

			if (MaxDelta < 0.001) {
				cout << endl << "END\n" << "N = " << n << ", " << " MaxDelta = " << MaxDelta << " < 10e3" << endl;
				break;
			}

			if (n % 10 == 0) {

				cout << "\tN = " << n << ", " << " MaxDelta = " << MaxDelta << endl;

				fout << MaxDelta << endl;
			}
	}
	
	return 0;

}