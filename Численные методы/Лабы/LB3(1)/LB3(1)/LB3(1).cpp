#include <iostream>
#include <vector>
#include <fstream>
const int n = 6;
using namespace std;
vector<double> gauss(vector<vector<double>> k, vector<double> R);

vector<double> gauss(vector <vector<double>> k, vector<double> R)
{
	double max;
	vector<double> x;
	x.resize(n);
	int c, index;
	const double eps = 0.00001;
	c = 0;
	while (c < n)
	{
		max = abs(k[c][c]);//поиск строки с максимальным k[c][c]
		index = c;
		for (int i = c + 1; i < n; i++)
		{
			if (abs(k[i][c]) > max)
			{
				max = abs(k[i][c]);
				index = i;
			}
		} 
		
		for (int j = 0; j < n; j++)//перестановка строк
		{
			double temp = k[c][j];
			k[c][j] = k[index][j];
			k[index][j] = temp;
		}
		double temp = R[c];
		R[c] = R[index];
		R[index] = temp;
		cout<< endl;
		cout << R[c]<<"\t"<<R[index] << endl;
			
		for (int i = c; i < n; i++)//нормализация уравнений
		{
			double temp = k[i][c];
			if (abs(temp) < eps) continue;//для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				k[i][j] = k[i][j] / temp;
			R[i] = R[i] / temp;
			if (i == c) continue;//уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)//вычитание строк
				k[i][j] = k[i][j] - k[c][j];
			R[i] = R[i] - R[c];
		}
		c++;
	}

	for (c = n - 1; c >= 0; c--)// обратная подстановка
	{
		x[c] = R[c];
		for (int i = 0; i < c; i++)
			R[i] = R[i] - k[i][c] * x[c];
	}

	return x;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	vector<double> C;
	vector<double> x{ -2.25 , -1.65, -1.05, -0.45, 0.15, 0.75 };
	vector<double> y{ 0.29422, 0.18737, -1.0215, -5.4471, -1.4440, 2.5873 };
	vector<double> R{ 1,1,1,1,1,1 };
	vector<vector<double>> A;
	A.resize(n);

	cout << "Задание 1 " << endl;
	cout << "Исходная матрица : " << endl;
	for (int i = 0; i < n; i++)//Вывод системы уравнений
	{
		A[i] = { -x[i], -pow(x[i],2), y[i], x[i] * y[i], y[i] * pow(x[i], 2), y[i] * pow(x[i], 3) };
		for (int j = 0; j < n; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << " = " << R[i] << endl;
	}

	cout << "-------------------------------------------------------" << endl;

	C = gauss(A, R);

	for (int i = 0; i < n; i++)
	{
		cout << "C[" << i << "] = " << C[i] << endl;
	}

	ofstream fout;
	fout.open("Result.txt");
	for (int i = 0; i < n; i++)
	{
		fout << x[i] << "  " << y[i] << endl;
	}
	fout.close();

	ofstream out;
	out.open("Result1.txt");
	out << 1 << endl;
	for (int i = 0; i < n; i++)
	{
		out << C[i] << endl;
	}
	out.close();

	return 0;

}


