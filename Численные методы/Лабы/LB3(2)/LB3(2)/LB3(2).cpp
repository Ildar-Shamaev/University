#include <iostream>
#include <vector>
#include <fstream>
const int n = 6;
using namespace std;

vector<double> Gauss(vector <vector<double>> k, vector<double> R)
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

void Print(vector<vector<double>> A)//вывод матрицы
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << "-------------------------------------------------------" << endl;
}

vector<vector<double>> LU(vector<vector<double>> U)
{
	vector<vector<double>> L;
	L.resize(n);
	for (int i = 0; i < n; i++) {
		L[i].resize(n * 2 + 1);
	}

	for (int i = 0; i < n; i++)
	{
		double max = 0;

		max = U[i][i];//k=U
		for (int j = i; j < n; j++)
		{
			L[j][i] = U[j][i] / max;
		}

		for (int j = i + 1; j < n; j++)
		{
			for (int z = i; z < n; z++)
			{
				U[j][z] -= U[i][z] * L[j][i];
			}
		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			L[i][j + n] = U[i][j];
		}
	}


	return L;
}

void Result(vector<double> C, vector<double> x, vector<double> y)
{
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
}


int main()
{
	setlocale(LC_ALL, "Russian");
	vector<double> C; vector<double> Y;

	vector<double> x{ -2.25 , -1.65, -1.05, -0.45, 0.15, 0.75 };
	vector<double> y{ 0.29422, 0.18737, -1.0215, -5.4471, -1.4440, 2.5873 };
	vector<double> R{ 1,1,1,1,1,1 };

	vector<vector<double>> A; A.resize(n);
	vector<vector<double>> lu;
	cout << "Задание 2" << endl;
	cout << "Исходная матрица :" << endl;
	vector<vector<double>> L; vector<vector<double>> U;
	L.resize(n); U.resize(n);
	for (int i = 0; i < n; i++)
	{
		L[i].resize(n); U[i].resize(n);
	}


	for (int i = 0; i < n; i++)
	{
		A[i] = { -x[i], -pow(x[i],2), y[i], x[i] * y[i], y[i] * pow(x[i], 2), y[i] * pow(x[i], 3) };
		for (int j = 0; j < n; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << " = " << R[i] << endl;
	}

	cout << "-------------------------------------------------------" << endl;

	lu = LU(A);



	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			L[i][j] = lu[i][j];
		}

		for (int j = n; j < 2 * n; j++)
		{
			U[i][j - n] = lu[i][j];
		}
	}

	cout << "L - матрица :" << endl;
	Print(L);
	cout << "U - матрица :" << endl;
	Print(U);

	Y = Gauss(L, R);
	C = Gauss(U, Y);

	for (int i = 0; i < n; i++)
	{
		cout << "C[" << i << "] = " << C[i] << endl;
	}

	Result(C, x, y);

	return 0;

}
