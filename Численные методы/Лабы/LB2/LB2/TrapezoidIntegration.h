#pragma once
#include "IntegrationMethod.h"

class TrapezoidIntegration : public IntegrationMethod
{
public:
	TrapezoidIntegration(double a = 0, double b = 0, unsigned n = 0, double delta = 0,
		double realResult = 0, double(*functionPtr)(double) = [](double x) { return 0.0; });
	virtual double calculateMethod();
	double calculateMethod(unsigned n);

private:

};