#pragma once
#include "IntegrationMethod.h"

class MonteCarlo : public IntegrationMethod
{
public:
	MonteCarlo(double a, double b, unsigned n, double delta, double realResult, double(*functionPtr)(double));
	virtual double calculateMethod() override;

	void checkExpectationFromN(unsigned points, unsigned repetitions, std::string fileName);

	double expectedValue(unsigned k);

private:


};
