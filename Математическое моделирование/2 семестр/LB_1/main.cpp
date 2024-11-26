#define _USE_MATH_DEFINES
#include <iostream>
#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>

using namespace std;

#define S 50
#define N 1000

struct double_pair
{
    double first;
    double second;
};

struct double_int
{
    int value;
    double key;
};

double_int find(int num, double_int* mas, int max)
{
    double_int result;
    for (int i = 0; i < max; i++)
    {
        if (mas[i].value == num)
        {
            if (num != 3)
            {
                result.key = mas[i].key;
                result.value = num;
                return result;
            }
            else
            {
                if (mas[i].key > double(3.55))
                {
                    result.key = mas[i].key;
                    result.value = num;
                    return result;
                }
            }
        }
    }
}

int main()
{
    int SpecialSolutionCount[] = { 2,4,8,16,3 };
    double begin = 3, end = 4, x = double(2.0 / 3.0), eps = 0.0001, step = 0.000001;
    double_int * solution_num = new double_int[1000000];
    int solution_count = 0;
    int tmp = 0;
    double func = x;
    double i = begin;
    while (i < end)
    {
        int unique_count = 0;
        double_int* unique_solutions = new double_int[100];
        for (int j = N; j > 1; j--)
        {
            func = func * i * (1.0 - func);
            if (j < S)
            {
                bool unique = true;
                for (int k = 0; k < unique_count; k++)
                {
                    if (abs(unique_solutions[k].key - func) < eps)
                    {
                        unique = false;
                        unique_solutions[k].value++;
                    }
                }
                if (unique)
                {
                    unique_solutions[unique_count].key = func;
                    unique_solutions[unique_count].value = 1;
                    unique_count++;
                }
            }
        }
        solution_num[tmp].key = i;
        solution_num[tmp].value = unique_count;
        tmp++;
        delete[] unique_solutions;
        i += step;
    }
    ofstream out;
    out.open("bif.csv");
    for (int i = 0; i < tmp; i++)
    {
        out << solution_num[i].key << ";" << solution_num[i].value << endl;
    }
    out.close();

    double_int* first_result = new double_int[5];
    for (int i = 0; i < 5; i++)
    {
        first_result[i] = find(SpecialSolutionCount[i], solution_num, tmp);
    }
    for (int i = 0; i < 5; i++)
    {
        cout << first_result[i].value << " " << first_result[i].key << " ";
        if (i > 1 && i <= 3)
        {
            double diff = (first_result[i - 1].key - first_result[i - 2].key) / (first_result[i].key - first_result[i - 1].key);
            cout << diff;
        }
        cout << endl;
    }
    _getch();
    return 0;
}