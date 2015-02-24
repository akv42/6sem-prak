#define _USE_MATH_DEFINES
#define L (10.0)
//define NUM в fft.h


#include "puasson.h"
#include "fft.h"
#include "thomas.h"

#include <iostream>
#include <fstream> 
#include <cmath>

using namespace std;

double f(const double &x, const double &y)
{
	return (sin(M_PI * x) * sin(M_PI * y));
}

void puas()
{
	//gO.Ogle: Лабораторная работа Применение быстрого преобразования Фурье для решения задачи распространения тепла в квадратной пластине
	//WARNING: туча очепяток в той работе <---
	cout.precision(5);
	double step = L / NUM;

	double **psson = new double*[NUM + 1];
	for (int i = 0; i < NUM; i++) {
		psson[i] = new double[NUM];
	}
	double *th_base = new double[NUM];
	double *ps_base = new double[NUM];
	double *base_f = new double[NUM];
	double *th_const = new double[NUM];
	double *th_c = new double[NUM];
	
	for (int i = 0; i < NUM; i++) {
		th_const[i] = 1;
	}

	double x = 0.0;
	double y = 0.0;
	for (int i_x = 0; i_x < NUM + 1; i_x++) {
		//БПФ для коэффициентов
		for (int i_y = 0; i_y < NUM; i_y++) {
			base_f[i_y] = f(x, y) * step * sqrt(2.0 / L);
			y += step;
		}
		th_base = fft_man(base_f);
		//подготовка к прогонке
		for (int i = 0; i < NUM; i++) {
			th_base[i] *= step * step;
			th_c[i] = 2.0 + 4.0 / step / step * pow(sin(M_PI * i_x * step / 2.0), 2);
		}
		//прогонка
		solveMatrix(NUM - 1, th_const, th_c, th_const, th_base, ps_base, 0.0, 0.0, 0.0, 0.0);
		//БПФ для решений
		for (int i_y = 0; i_y < NUM; i_y++) {
			ps_base[i_y] *= sqrt(2.0 / L);
		}
		psson[i_x] = fft_man(ps_base);
		for (int i_y = 0; i_y < NUM; i_y++) {
			cout << fixed << psson[i_x][i_y] << ' ';
		}
		x += step;
		y = 0.0;

		cout << fixed << y << endl; //не считаем границу (она все равно 0, а у меня fft поляжет)
	}

	for (int i = 0; i < NUM; i++) {
		delete[] psson[i];
	}
	delete[] psson;
	delete[] th_base;
	delete[] ps_base;
	delete[] base_f;
	delete[] th_const;
	delete[] th_c;
}