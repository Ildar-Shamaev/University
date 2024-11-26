#pragma once
#include "IntegrationMethod.h"
#include "TrapezoidIntegration.h"
#include <vector>

using namespace std;

class OptimalNodes : public IntegrationMethod
{
public:
	OptimalNodes(double a, double b, unsigned n, unsigned q, double realResult, double delta, double(*functionPtr)(double),
		double(*secondDerivative)(double), vector<double>& zerosThirdDerivative);

	virtual double calculateMethod() override;

private:
	double maxFunctionOnSegment(double x_i, double x_j);
	double convertSegmentToGivenInterval(double point);
	double convertSegmentForMethod(double point);
	double(*secondDerivative)(double);
	vector<double> zerosThirdDerivative;
	vector<pair<double, double>> b_l;
	vector<double> A_l;
	unsigned q;

	TrapezoidIntegration trapezoid;
};
