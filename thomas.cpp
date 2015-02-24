#include "thomas.h"

#include <iostream>

using namespace std;

void solveMatrix(int N, double *a, double *c, double *b, double *f, double *x, double k1, double k2, double m1, double m2)
{
	double *p = new double[N + 1];
	double *q = new double[N + 1];
	p[1] = k1;
	q[1] = m1;
	for (int i = 1; i < N; i++)
	{
		p[i + 1] = b[i] / (c[i] - a[i] * p[i]);
		q[i + 1] = (a[i] * q[i] + f[i]) / (c[i] - a[i] * p[i]);
	}

	x[N] = (m2 + k2 * q[N]) / (1 - k2 * p[N]);

	for (int i = N; i > 0; i--) {
		x[i - 1] = p[i] * x[i] + q[i];
	}

	delete[] p;
	delete[] q;
}