#pragma once
#include <string>

using namespace std;

class IntegrationMethod
{
public:
	IntegrationMethod(double a, double b, unsigned n, double realResult, double delta, double(*functionPtr)(double));

	unsigned minError(string fileName);

	inline void setAB(double a, double b) { this->a = a; this->b = b; }

	virtual double calculateMethod() { return 0.0; };
	double calculateMethod(unsigned n);

	void setGrid(unsigned n);




protected:
	string firstSchedule;
	double a;
	double b;
	unsigned n;
	double realResult;
	double (*functionPtr)(double);

	double delta;
};