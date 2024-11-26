#include "MonteCarlo.h"
#include <time.h>
#include <random>	
#include <fstream>

using namespace std;

MonteCarlo::MonteCarlo(double a, double b, unsigned n, double delta, double realResult, double(*functionPtr)(double)) : IntegrationMethod(a, b, n, delta, realResult, functionPtr)
{
	srand(time(0));
}

double MonteCarlo::calculateMethod()
{
	//������� �������� �����-����� ��������������
	double u = 0.;//��������� �������� �������������� �� ������� ��������������
	double sum = 0.;
	double x = 0.;
	double value = 0.;
	for (int i = 0; i < n; i++) {
		u = double(rand()) / double(RAND_MAX);//��������� ��������
		x = a + (b - a) * u;
		sum += functionPtr(x);
		value = sum / n;//���������� ������� ��� ���������� ������ ���������
	}
	return ((b - a) * value);
}

void MonteCarlo::checkExpectationFromN(unsigned points, unsigned repetitions, string fileName)//���.�������� �� ���������� N(����� ������ �������� ������) 
{
	ofstream outFile(fileName);

	for (unsigned i = 1; i <= points; ++i) {
		setGrid(i);
		outFile << i << " " << expectedValue(repetitions) << "\n";
	}
	outFile.close();

}


double MonteCarlo::expectedValue(unsigned k)//��������� ��������
{
	double expectRes = 0;
	for (int i = 0; i < k; i++)
	{
		expectRes += fabs(realResult - calculateMethod());
	}

	return expectRes / k;

}
