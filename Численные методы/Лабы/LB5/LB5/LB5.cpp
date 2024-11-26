#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
using namespace std;
void Out(vector<vector<double>> vec)//Вывод матрици
{
	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec[i].size(); j++)
		{
			cout << vec[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
void Out(vector<double> vec)//вывод вектора
{
	for (int i = 0; i < vec.size(); i++)
	{
		cout << vec[i] << "\n";
	}
	cout << endl;
}

vector<double> Zap(int N)//Заполнить вектор начального значения
{
	vector<double> vec; 
	vec.resize(N);
	for (int i = 0; i < N; i++)
	{
		vec[i] = 0.5;//Почему 0.5!!!!!!!!!!!!!!!!!!
	}
	return vec;
}

double Norm(vector<double> vec, vector<double> vec1)//в цикле 2 вектора минус и его норма
{
	double norm = 0;
	for (int i = 0; i < vec.size(); i++)
	{

		norm += pow((vec[i] - vec1[i]), 2);
	}
	return sqrt(norm);
}

vector<vector<double>> Lent(int N, int l)//Матрици 3 ленточные
{
	vector<double> q;
	q = { 1.1,2,10 };
	vector<vector<double>> vec;
	vec.resize(N);
	for (int i = 0; i < N; i++)
	{
		vec[i].resize(3 * N + 1);
	}
		

	for (int i = 0; i < N; i++)//заполнение побочных диаг
	{
		int left = i - l;
		if (left < 0)
		{
			left = 0;
		}
		int right = i + l + 1;
		if (right > N)
		{
			right = N;
		}
		for (int j = left; j < right; j++)//Заполнение диагонали случайными числами
		{
			if (i == j) continue;
			vec[i][j] = (rand() % 1001 - 50) / 1000.0;
		}
	}
	//Out(vec);
	for (int z = 0; z < q.size(); z++)//заполнение диагоналейй с q
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = z * N; j < z * N + N; j++)
			{
				if ((i + z * N) == j)
				{
					double Sum = 0;
					for (int shg = 0; shg < N; shg++)
					{
						if (shg == i) continue;
						Sum += abs(vec[i][shg]);
					}
					vec[i][j] = q[z] * Sum;
					continue;
				}
				vec[i][j] = vec[i][j - z * N];
			}
		}
	}
	//Out(vec);
	for (int i = 0; i < N; i++)//вектор точного решения СЛАУ
	{
		vec[i][N * 3] = (rand() % 1001 - 50) / 1000.0;
	}
	//Out(vec);
	return vec;
}

vector<vector<double>> Vprav(vector<vector<double>> vec, int N)//Вектора правые для соответсвующей ленточной матрици
{
	vector<vector<double>> b;
	b.resize(N);
	for (int i = 0; i < N; i++) b[i].resize(3);

	for (int z = 0; z < 3; z++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = z * N; j < z * N + N; j++)
			{
				b[i][z] += vec[i][j] * vec[j - z * N][3 * N];//вычисление 3 векторов правой части системы
			}
		}
	}

	return b;
}

vector<vector<double>> Transp(vector<vector<double>> vec, int N)//транспонирование 3 матриц
{
	for (int i = 0; i < N; i++)
	{
		for (int j = i; j < N; j++)
		{
			if (i == j) continue;
			double per = vec[i][j];
			vec[i][j] = vec[j][i];
			vec[i][j + N] = vec[j][i];
			vec[i][j + 2 * N] = vec[j][i];
			vec[j][i] = per;
			vec[j][i + N] = per;
			vec[j][i + 2 * N] = per;
		}
	}
	//Out(vec);
	return vec;
}

vector<vector<double>> Proizved(vector<vector<double>> vec, vector<vector<double>> vect, int N)//произведение матриц для придания симметричности
{
	vector<vector<double>> Rezvec;
	Rezvec.resize(N);
	for (int i = 0; i < N; i++)
	{
		Rezvec[i].resize(3 * N);
	}
	for (int z = 0; z < 3; z++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = z * N; j < z * N + N; j++)//По координатам результата
			{
				double Sum = 0;
				for (int jy = z * N; jy < z * N + N; jy++)//сложение элементов матрицы
				{
					Sum += vect[i][jy] * vec[jy - z * N][j];
				}
				Rezvec[i][j] = Sum;
			}
		}
	}
	//Out(Rezvec);
	return Rezvec;
}

vector<vector<double>> proizvedB(vector<vector<double>> vect, vector<vector<double>> b, int N)//соответствующие правые вектора для соответствующих симметричых матриц
{
	vector<vector<double>> B;
	B.resize(N);
	for (int i = 0; i < N; i++)
	{
		B[i].resize(3);
	}
	for (int z = 0; z < 3; z++)
	{
		for (int i = 0; i < N; i++)
		{
			double Sum = 0;
			for (int j = z * N; j < z * N + N; j++)
			{
				Sum += vect[i][j] * b[j - z * N][z];
			}
			B[i][z] = Sum;
		}
	}
	//Out(B);
	return B;

}

int Icobi(vector<vector<double>> vec, vector<double> nach, int N, double eps)//Метод Якоби
{
	int it = 0;

	//Преобразование B

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j) continue;
			vec[i][j] /= ((-1) * vec[i][i]);
		}
		vec[i][N] /= vec[i][i];
		vec[i][i] = 0;
	}

	vector<double> rez; 
	rez = nach;
	do
	{
		nach = rez;
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			for (int j = 0; j < N; j++)
			{
				sum += vec[i][j] * nach[j];

			}
			rez[i] = sum + vec[i][N];
		}
		it++;
	} while (Norm(nach, rez) > eps);//среднеквадратичное отклонение
	//Out(rez);

	ofstream Outv;
	Outv.open("Vec.txt");
	for (int i = 0; i < N; i++)
	{
		Outv << rez[i] << " " << endl;
	}
	return it;
}

int SOR(vector<vector<double>> vec, vector<double> nach, int N, double eps, double w)
{
	int it = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j) continue;
			vec[i][j] /= ((-1) * vec[i][i]);
		}
		vec[i][N] /= vec[i][i];
		vec[i][i] = 0;
	}

	vector<double> rez; rez = nach;
	do
	{
		rez = nach;
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			for (int j = 0; j < N; j++)
			{
				sum += vec[i][j] * nach[j];
			}
			nach[i] = sum + vec[i][N];
		}

		for (int i = 0; i < N; i++)
		{
			nach[i] = nach[i] * w + (1 - w) * rez[i];//SOR
		}
		it++;
	} while (Norm(nach, rez) > eps);
	//if (abs(w - 1) < 0.01) { cout << endl; cout << "Вектор при w = 1" << endl; Out(nach); }
	return it;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	srand(time(0));
	int N = 75, l = 5; double eps = 0.000001;
	vector<vector<double>>vec, vect, vecRez;
	vector<vector<double>>b;
	vector<double> q;
	q = { 1.1,2,10 };

	vec = Lent(N, l);
	//Out(vec);

	b = Vprav(vec, N);
	
	vect = Transp(vec, N);

	vecRez = Proizved(vec, vect, N);
	for (int i = 0; i < N; i++)
	{
		vecRez[i].push_back(vec[i][3 * N]);
	}

	b = proizvedB(vect, b, N);//b*
	//Вывод 3 матриц
	
	vector<vector<double>> Rez; 
	Rez.resize(N);
	for (int i = 0; i < N; i++)
	{
		Rez[i].resize(N + 1);
	}
	/*
	for (int z = 0; z < 3; z++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = z * N; j < N + z * N; j++)
			{
				Rez[i][j - z * N] = vecRez[i][j];
			}
			Rez[i][N] = b[i][z];

		}
		//cout << q[z] << endl;
		//Out(Rez);
	}
	*/
	cout << "Метод Якоби" << endl;

	//Задание 2

	vector<double> nach;
	nach = Zap(N);
	//cout << "Начальный вектор: " << endl; Out(nach);
	ofstream Outt;
	Outt.open("Matrix.txt");
	ofstream Outv;
	Outv.open("Vec.txt");

	for (int z = 0; z < 3; z++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = z * N; j < N + z * N; j++)
			{
				Rez[i][j - z * N] = vecRez[i][j];
				Outt << Rez[i][j - z * N] << " ";
			}
			Rez[i][N] = b[i][z];
			Outt << Rez[i][N] << endl;
		}
		Outt << endl;
		cout << "q= " << q[z] << endl;
		cout << "Количнство итераций = " << Icobi(Rez, nach, N, eps) << endl;

	}
	Outt.close();


	////Задание 3
	cout << "Метод SOR" << endl;

	for (int z = 0; z < 3; z++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = z * N; j < N + z * N; j++)
			{
				Rez[i][j - z * N] = vecRez[i][j];
			}
			Rez[i][N] = b[i][z];
		}
		cout << "q= " << q[z] << endl;

		cout << "W		Количество итераций" << endl;

		for (double w = 0.1; w < 2; w += 0.1)
		{
			cout << w << "		" << SOR(Rez, nach, N, eps, w) << endl;
		}
	}
}