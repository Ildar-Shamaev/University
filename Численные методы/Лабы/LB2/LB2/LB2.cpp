#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "MonteCarlo.h"
#include "OptimalNodes.h"
#include "TrapezoidIntegration.h"

#define realResult 0.612889141
#define delta 0.001
#define a 0.
#define b 2
#define n 100

using namespace std;

double Function(double x);
double FunctionSecondDerivative(double x);

int main() {
	
	setlocale(0, "Russian");
	
	MonteCarlo montecarlo{ a, b, n, realResult, delta, Function };
	cout << setprecision(10) << "Интеграл по Монтекарло = " << montecarlo.calculateMethod() << "\n";
	cout << setprecision(10) << "Математическое ожидание Монтекарло = " << montecarlo.expectedValue(100) << "\n";
	montecarlo.checkExpectationFromN(1000, 100, "checkmateWaiting.dat");
	
	cout << "\n\n";

	vector<double> ZerosThirdDerivative{ 0.37918, 1.84659 };

	OptimalNodes optimalTrapezoid{ a, b, n, 10, realResult, delta, Function, FunctionSecondDerivative, ZerosThirdDerivative };
	cout << setprecision(10) << "Метод оптимального распределения узлов трапеции = " << optimalTrapezoid.calculateMethod() << "\tОшибка: " << fabs(optimalTrapezoid.calculateMethod() - realResult) << "\n";
	cout << "Минимальное n метода оптимального распределения узлов трапеции = " << optimalTrapezoid.minError("optimalTrapezoidError.dat") << "\n";
	
}

double Function(double x) {
	return atan(x)/(1+x*x);
}
double FunctionSecondDerivative(double x) {

	return (-6*x-2*atan(x)+6*atan(x)*x*x)/pow(1+x*x,3);
}


