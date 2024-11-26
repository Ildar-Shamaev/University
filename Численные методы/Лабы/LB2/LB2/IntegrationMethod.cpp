#include "IntegrationMethod.h"
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

IntegrationMethod::IntegrationMethod(double a, double b, unsigned n, double realResult, double delta, double(*functionPtr)(double))
{
	this->a = a;
	this->b = b;
	this->n = n;
	this->realResult = realResult;
	this->delta = delta;
	this->functionPtr = functionPtr;
}

unsigned IntegrationMethod::minError(string fileName)
{
	ofstream out(fileName);
	double error = 0.0;
	n = 0;
	do
	{
		n++;
		double temp = calculateMethod();
		//error = fabs(realResult - calculateMethod());
		double result = fabs(realResult - temp);
		out << n << " " << result << "\n";
		error = result;
		
	} while (error > delta);
	out.close();
	return n;
}
/*
double IntegrationMethod::calculateMethod(unsigned n)
{
	unsigned lastN = this->n;
	setGrid(n);

	double result = calculateMethod();
	setGrid(lastN);

	return result;
}
*/
void IntegrationMethod::setGrid(unsigned n)
{
	this->n = n;

}



