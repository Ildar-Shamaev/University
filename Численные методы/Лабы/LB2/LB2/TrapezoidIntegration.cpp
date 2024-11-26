#include "TrapezoidIntegration.h"
#include <math.h>
#include <iostream>

TrapezoidIntegration::TrapezoidIntegration(double a, double b, unsigned n, double delta, double realResult,
	double(*functionPtr)(double)) : IntegrationMethod(a, b, n, delta, realResult, functionPtr)
{
}

double TrapezoidIntegration::calculateMethod()
{
	double sum = 0.0;
	double h = (b - a) / n;
	double x_i = 0.0;
	double x_j = 0.0;

	for (unsigned i = 0; i < n; i++)
	{
		x_i = a + h * i;
		x_j = a + h * (i + 1);
		sum += (functionPtr(x_i) + functionPtr(x_j));
	}
	sum = sum * 0.5 * h;//формула трапеции
	return sum;
}

double TrapezoidIntegration::calculateMethod(unsigned n)
{
	setGrid(n);
	return calculateMethod();
}
