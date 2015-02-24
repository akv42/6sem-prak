#pragma comment (lib,"libfftw3-3.lib")

#define _USE_MATH_DEFINES

#include "fft.h"

#include <fstream> 
#include <iostream>
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <cstdio>
#include <ctime>

using namespace std;

ofstream fft_out;

void file_out(double *result)
{
	for (int i = 0; i < NUM; i++) {
		fft_out << result[i] << ' ';
	}
	fft_out << endl << endl;
}

int order_2()
{
	int ord = 0;
	int num = NUM;
	while (num != 1) {
		num /= 2;
		ord++;
	}
	return ord;
}

void transform(complex<double> *fft_data, int num_tr)
{
	complex<double> *tmp_0 = new complex<double>[num_tr / 2];
	complex<double> *tmp_1 = new complex<double>[num_tr / 2];

	for (int i = 0; i < num_tr; i += 2) {
		tmp_0[i / 2] = fft_data[i];
		tmp_1[i / 2] = fft_data[i + 1];
	}

	for (int i = 0; i < num_tr; i += 2) {
		fft_data[i / 2] = tmp_0[i / 2];
		fft_data[i / 2 + num_tr / 2] = tmp_1[i / 2];
	}

	//рекурсия
	if (num_tr / 4 != 1) {
		transform(fft_data, num_tr / 2);
		transform(fft_data + (num_tr / 2), num_tr / 2);
	}

	delete[] tmp_0;
	delete[] tmp_1;
}

void func(fftw_complex *fft_data)
{
	complex<double> m_i(0.0, 1.0);
	complex<double> cur(0.0, 0.0);
	int t = -NUM / 2;
	for (int i = 0; i < NUM; i++) {
		cur = sin(M_PI * t);
		fft_data[i][0] = real(cur);
		fft_data[i][1] = imag(cur);
		t++;
	}
}

void func(complex<double> *fft_base)
{
	complex<double> m_i(0.0, 1.0);
	int t = -NUM / 2;
	for (int i = 0; i < NUM; i++) {
		fft_base[i] = sin(t);
		t++;
	}
}

complex<double> fft_exp(int up, int down)
{
	complex<double> m_i(0.0, 1.0);
	return exp(-m_i * 2.0 * M_PI * (double)up / (double)down);
}

void fft_make(complex<double> *fft_data, int down, int num_tr)
{
	complex<double> *tmp = new complex<double>[NUM];

	for (int i = 0; i < num_tr; i += 2) {
		for (int up = 0; up < down / 2; up++) {
			tmp[up + i * down / 2] = fft_data[up + i * down / 2] + fft_exp(up, down) * fft_data[up + (i + 1) * down / 2];
		}
		for (int up = 0; up < down / 2; up++) {
			tmp[up + (i + 1) * down / 2] = fft_data[up + i * down / 2] - fft_exp(up, down) * fft_data[up + (i + 1) * down / 2];
		}
	}

	for (int i = 0; i < NUM; i++) {
		fft_data[i] = tmp[i];
	}
	
	if (down != NUM) {
		fft_make(fft_data, 2 * down, num_tr / 2);
	}

	delete[] tmp;
}

void fft_man()
{
	//подробно: http://www.dsplib.ru/content/thintime/thintime.html
	int order = order_2();
	
	complex<double> *fft_data = new complex<double>[NUM];

	//массивы результатов
	double *result = new double[NUM];

	//заполняем массив значениями функции
	func(fft_data);

	//рекурсивная перестановка на "четный-нечетный"
	transform(fft_data, NUM);

	//рекурсивное формирование спектра бабочкой
	fft_make(fft_data, 2, NUM);

	//строим функцию полученную (как модуль комплексного числа)
	for (int i = 0; i < NUM; i++) {
		result[i] = abs(fft_data[i]);
	}

	file_out(result);

	delete[] result;
	delete[] fft_data;
}

double* fft_man(double *base)
{
	//подробно: http://www.dsplib.ru/content/thintime/thintime.html
	int order = order_2();

	complex<double> *fft_data = new complex<double>[NUM];

	//массивы результатов
	double *result = new double[NUM];

	//заполняем массив значениями функции
	for (int i = 0; i < NUM; i++) {
		fft_data[i] = base[i];
	}

	//рекурсивная перестановка на "четный-нечетный"
	transform(fft_data, NUM);

	//рекурсивное формирование спектра бабочкой
	fft_make(fft_data, 2, NUM);

	//берем мнимую часть (синусы)
	for (int i = 0; i < NUM; i++) {
		result[i] = imag(fft_data[i]);
	}

	delete[] fft_data;

	return result;
}

void fft_auto()
{
	fftw_complex *fft_res = new fftw_complex[NUM];
	fftw_complex *fft_base = new fftw_complex[NUM];

	//массивы результатов
	double **result = new double*[2];
	result[0] = new double[NUM];
	result[1] = new double[NUM];
	
	//заполняем массив значениями функции
	func(fft_base);

	//строим функцию исходную (как модуль комплексного числа)
	for (int i = 0; i < NUM; i++) {
		result[0][i] = sqrt(fft_base[i][0] * fft_base[i][0] + fft_base[i][1] * fft_base[i][1]);
	}

	//вычисляем спектр
	fftw_plan plan = fftw_plan_dft_1d(NUM, fft_base, fft_res, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	//строим функцию полученную (как модуль комплексного числа)
	for (int i = 0; i < NUM; i++) {
		result[1][i] = sqrt(fft_res[i][0] * fft_res[i][0] + fft_res[i][1] * fft_res[i][1]);
	}

	fftw_destroy_plan(plan);

	file_out(result[1]);

	delete[] result[0];
	delete[] result[1];
	delete[] result;

	delete[] fft_res;
	delete[] fft_base;
}

void fft_test()
{
	fft_out.open("fft_out.txt");
	clock_t t_begin, t_end;
	t_begin = clock();
	fft_man();
	t_end = clock();
	cout << "fft_man: " << (double)(t_end - t_begin) / CLOCKS_PER_SEC << " sec" << endl;
	t_begin = clock();
	fft_auto();
	t_end = clock();
	cout << "fft_auto: " << (double)(t_end - t_begin) / CLOCKS_PER_SEC << " sec" << endl;
	fft_out.close();
}