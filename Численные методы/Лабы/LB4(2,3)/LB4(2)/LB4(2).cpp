#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
//#include<iomanip>

using namespace std;
const double eps = 10e-6;
//2 ЗАДАНИЕ
double F_x(double x, double y)
{
	return -cos(x) * cos(y) / 6.0;
}
double F_y(double x, double y) {

	return 1 + sin(x) * sin(y) * 0.5;
}

//3 ЗАДАНИЕ
double F1(double x, double y)
{

	return 6 * x + cos(y) * cos(x);
}
double F2(double x, double y)
{

	return  2 * y - 2 - sin(x) * sin(y);
}
double F1dx(double x, double y)
{
	return 6 - cos(y) * sin(x);
}
double F1dy(double x, double y)
{
	return -cos(x) * sin(y);
}
double F2dx(double x, double y)
{

	return -sin(y) * cos(x);
}

double F2dy(double x, double y)
{
	return 2 - sin(x) * cos(y);
}

double det(double x, double y)//определитель
{
	double det = F1dx(x, y) * F2dy(x, y) - F1dy(x, y) * F2dx(x, y);
	return det;
}


double FuncXY(double x, double y){//заданная функция

	return x * x + y * y + (x + y + 1) * tan(x + y);
}

void SimpleIterationsMethod(double x0, double y0)
{
	int iteration = 0;
	double x_k1 = x0;
	double x_k2 = x_k1;
	double y_k1 = y0;
	double y_k2 = y_k1;
	double delta = 1;
	while (delta >= eps)
	{
		iteration++;
		x_k2 = x_k1;
		y_k2 = y_k1;
		x_k1 = F_x(x_k2, y_k2);
		y_k1 = F_y(x_k2, y_k2);
		delta = fmax(abs(x_k2 - x_k1), abs(y_k2 - y_k1));
		if (delta >= eps) {
			cout << iteration << ": F(" << x_k1 << "; " << y_k1 << ") = " << FuncXY(x_k1, y_k1) << "\tdelta: " << delta << " > " << 0.000001 << endl;
		}
		else
			cout << iteration << ": F(" << x_k1 << "; " << y_k1 << ") = " << FuncXY(x_k1, y_k1) << "\tdelta: " << delta << " < " << 0.000001 << endl;
	}
	cout << endl;
	cout << "Точка минимума: (" << x_k1 << "; " << y_k1 << ")" << endl;
	cout << " F(" << x_k1 << "; " << y_k1 << ") = " << FuncXY(x_k1, y_k1) << endl;
	cout << " Количество итераций: " << iteration << endl;

}


void NewtonMethodExact(double x0, double y0)
{
	vector<double> root;
	int iteration = 0;
	double x_k = x0;
	double x_k1 = x_k;
	double y_k = y0;
	double y_k1 = y_k;
	double delta = 1;
	//Метод основан на матрице обратной к матрице Якоби
	//Значение Якобиана
	while (delta >= eps)
	{
		x_k1 = x_k;
		y_k1 = y_k;

		x_k = x_k1 - (1.0 / det(x_k1, y_k1)) * (F2dy(x_k1, y_k1) * F1(x_k1, y_k1) - F1dy(x_k1, y_k1) * F2(x_k1, y_k1));
		y_k = y_k1 - (1.0 / det(x_k1, y_k1)) * (F1dx(x_k1, y_k1) * F2(x_k1, y_k1) - F2dx(x_k1, y_k1) * F1(x_k1, y_k1));

		iteration++;
		delta = fmax(abs(x_k1 - x_k), abs(y_k1 - y_k));
		if (delta >= eps) {
			cout << iteration << ": F(" << x_k << "; " << y_k << ") = " << FuncXY(x_k, y_k) <<"\tdelta: " << delta << " > " << 0.000001 << endl;
		}
		else
			cout << iteration << ": F(" << x_k << "; " << y_k << ") = " << FuncXY(x_k, y_k) <<"\tdelta: " << delta << " < " << 0.000001 << endl;
	}
	cout << endl;
	cout << "Точка минимума: (" << x_k1 << "; " << y_k1 << ")" << endl;
	cout << " F(" << x_k1 << "; " << y_k1 << ") = " << FuncXY(x_k1, y_k1) << endl;
	cout << " Количество итераций: " << iteration << endl;
}

void NewtonMethodApprox(double x0, double y0)
{
	vector<double> root;
	int iteration = 0;
	double x_k = x0;
	double x_k1 = x_k;
	double y_k = y0;
	double y_k1 = y_k;
	double delta = 1;
	double h = 0.0001;
	//Метод основан на матрице обратной к матрице Якоби
	//Значение Якобиана
	while (delta >= eps)
	{
		x_k1 = x_k;
		y_k1 = y_k;
		/*
		Основная сложность метода Ньютона заключается в обращении матрицы Якоби.
		Вводя обозначение Δx(k) = x(k + 1)−x(k) получаем СЛАУ для вычисления
		Δx(k) : ∂F(x(k))/∂x = −F(x(k))
		Тогда x(k + 1) = x(k) + Δx(k).
		*/
		double F1dx = (F1(x_k + h, y_k) - F1(x_k, y_k)) / h;
		double F1dy = (F1(x_k, y_k + h) - F1(x_k, y_k)) / h;
		double F2dx = (F2(x_k + h, y_k) - F2(x_k, y_k)) / h;
		double F2dy = (F2(x_k, y_k + h) - F2(x_k, y_k)) / h;
		double det = F1dx * F2dy - F1dy * F2dx;

		x_k = x_k1 - (1.0 / det) * (F2dy * F1(x_k1, y_k1) - F1dy * F2(x_k1, y_k1));
		y_k = y_k1 - (1.0 / det) * (F1dx * F2(x_k1, y_k1) - F2dx * F1(x_k1, y_k1));

		iteration++;
		delta = fmax(abs(x_k1 - x_k), abs(y_k1 - y_k));
		if (delta >= eps) {
			cout << iteration << ": F(" << x_k << "; " << y_k << ") = " << FuncXY(x_k, y_k) << "\tdelta: " << delta << " > " << 0.000001 << endl;
		}
		else
			cout << iteration << ": F(" << x_k << "; " << y_k << ") = " << FuncXY(x_k, y_k) << "\tdelta: " << delta << " < " << 0.000001 << endl;
	}
	cout << endl;
	cout << "Точка минимума: (" << x_k1 << "; " << y_k1 << ")" << endl;
	cout << " F(" << x_k1 << "; " << y_k1 << ") = " << FuncXY(x_k1, y_k1) << endl;
	cout << " Количество итераций: " << iteration << endl;
}
int main()
{
	setlocale(LC_ALL, "Russian");
	vector<double> root;
	double x0 = -1/2;
	double y0 = 0;

	SimpleIterationsMethod(x0, y0);
	cout << "\n____________________________\n";
	NewtonMethodExact(x0, y0);
	cout << "\n____________________________\n";
	NewtonMethodApprox(x0, y0);
	return 0;
}

