#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#define pi 3.14159
#define _USE_MATH_DEFINES

const double a = 0, b = 2;

using namespace std;

double f(double x)
{
	return atan(x) / (1 + x * x);
}
void Progonka(double* y, double* C, int n, double h)//4N неизвестных
{
	double* alpha = new double[n + 1];//переменная для предположения
	double* beta = new double[n + 1];//переменная для предположения
	double* F = new double[n + 1];//Функция для прогонки
	double z = 0;//xi предположения
	double c = 4 * h;
	double b = h;
	double a = h;
	int N = n - 1;

	for (int i = 1; i <= N; i++)
	{
		F[i] = (3 / h) * (y[i + 1] - 2 * y[i] + y[i - 1]);//уравнение для С (28)
	}

	alpha[1] = -c / b;
	beta[1] = F[1] / b;

	for (int i = 2; i <= N - 1; i++)
	{
		z = b + a * alpha[i - 1];//шаг
		alpha[i] = -c / z;
		beta[i] = (F[i] - a * beta[i - 1]) / z;
	}

	beta[N] = (F[N] - a * beta[N - 1]) / (b + a * alpha[N - 1]);
	C[n] = 0;
	C[N] = beta[N];

	for (int i = N - 1; i >= 1; i--)
	{
		C[i] = alpha[i] * C[i + 1] + beta[i];
	}

	C[0] = 0;
}

void coefPn(double* y, double* B, double* C, double* D, int n, double h)
{
	for (int i = 1; i <= n - 1; i++)
	{
		D[i] = (C[i] - C[i - 1]) / (3 * h);//(26)
	}

	D[n] = -C[n - 1] / (3 * h);//(26)

	for (int i = 1; i <= n - 1; i++)
	{
		B[i] = ((y[i] - y[i - 1]) / h) - (h / 3) * (C[i] + 2 * C[i - 1]);//(27)
	}

	B[n] = ((y[n] - y[n - 1]) / h) - (2 * h / 3) * C[n - 1];//(27)
}
double Pn(double xi, double* y, double* x, double* b, double* c, double* d, int n)//Сплайн-функция
{
	for (int i = 0; i <= n; i++)
	{
		if (xi >= x[i] && xi <= x[i + 1])
			return (y[i] + b[i + 1] * (xi - x[i]) + c[i] * pow((xi - x[i]), 2) + d[i + 1] *
				pow((xi - x[i]), 3));
	}
}
int main()
{
	int n;
	for (int n = 1; n < 49; n++)
	{
		double h = (b - a) / n;//сетка
		double* x = new double[n + 1];
		double* y = new double[n + 1];

		for (int i = 0; i <= n; i++)
		{
			x[i] = i * h;
			y[i] = f(x[i]);
		}

		double* B = new double[n + 1];
		double* C = new double[n + 1];
		double* D = new double[n + 1];
		Progonka(y, C, n, h);
		coefPn(y, B, C, D, n, h);
		double delta = 0.0;
		double maxDelta = 0.0;

		for (double i = 0; i <= b; i += h)
		{
			delta = abs(f(i) - Pn(i, y, x, B, C, D, n));//Вычисление погрешности
			if (maxDelta <= delta)
				maxDelta = delta;
		}

		cout << "delta[" << n << "] = " << maxDelta;
		if (maxDelta <= 0.001)
		{
			cout << " - Podhodit" << endl;
		}
	}
	return 0;
}
