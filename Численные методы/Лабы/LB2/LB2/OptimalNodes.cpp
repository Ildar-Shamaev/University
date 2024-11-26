#include "OptimalNodes.h"
#include <iostream>

using namespace std;

OptimalNodes::OptimalNodes(double a, double b, unsigned n, unsigned q, double realResult, double delta, double(*functionPtr)(double), double(*secondDerivative)(double),
	vector<double>& zerosThirdDerivative) : IntegrationMethod(a, b, n, realResult, delta, functionPtr)
{
	trapezoid = { a, b, n, delta, realResult, functionPtr };
	this->secondDerivative = secondDerivative;
	this->q = q;

	for (int i = 0; i < zerosThirdDerivative.size(); i++)
	{
		this->zerosThirdDerivative.push_back(convertSegmentForMethod(zerosThirdDerivative[i]));
	}
}

double OptimalNodes::calculateMethod()
{
	double h = 1. / q;//шаг

	double x_i = 0.;
	double x_j = 0.;
	double sum = 0.;
	double six_lambda = 0.;//оценка погрешности по элементраному отрезку интегрирования bl/Nl
	double answer = 0.;
	double checkN;
	double tempRes;


	for (int i = 0; i < q; ++i) {
		x_i = h * i;
		x_j = h * (i + 1);
		b_l.push_back(make_pair(x_i, x_j));//создаем пару(интервал)
		A_l.push_back(maxFunctionOnSegment(x_i, x_j));
		

		sum += (x_j - x_i) * pow(A_l.back(), 1. / 3.);
	}
	six_lambda = sum / n;
	double tPow;
	for (int i = 0; i < q; ++i) {
		trapezoid.setAB(convertSegmentToGivenInterval(b_l[i].first), convertSegmentToGivenInterval(b_l[i].second));

		tPow = pow(A_l[i], 1. / 3.);

		tempRes = (int)(((b_l[i].second - b_l[i].first)) * tPow / six_lambda);

		checkN = tempRes != 0 ? tempRes : 1;


		trapezoid.setGrid(checkN);

		answer += trapezoid.calculateMethod();

	}

	return answer;

}

double OptimalNodes::maxFunctionOnSegment(double x_i, double x_j)
{
	double functionZero;
	double functionA = fabs(secondDerivative(convertSegmentToGivenInterval(x_i)));
	double functionB = fabs(secondDerivative(convertSegmentToGivenInterval(x_j)));

	for (int i = 0; i < zerosThirdDerivative.size(); ++i) {

		if (x_i <= zerosThirdDerivative[i] && x_j >= zerosThirdDerivative[i]) //находит интервал ф, где 3 производная обращается в 0
		{
			functionZero = fabs(secondDerivative(convertSegmentToGivenInterval(zerosThirdDerivative[i])));
			return ((functionZero >= functionA && functionZero >= functionB) ? 
				functionZero : ((functionA >= functionB) ? functionA : functionB));
		}
	}
	return ((functionA >= functionB) ? functionA : functionB);
}

double OptimalNodes::convertSegmentToGivenInterval(double point)
{
	return ((b - a) * point + a);
}

double OptimalNodes::convertSegmentForMethod(double point)
{
	return ((point - a) / (b - a));
}
