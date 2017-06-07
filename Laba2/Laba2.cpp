// Laba2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <thread>
#include <fstream>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <ctime> 

using namespace std;

void Func1(double *vec1, double *vec2, double *scalar_mult, int idx, int n)
{
	double s3 = 0;
	int k;
	for (k = idx; k < n; k = k + 2)
	{
		s3 = s3 + vec1[k] * vec2[k];
	}
	scalar_mult[0] = scalar_mult[0] + s3;
}
void Func(double **A, double *S, double *Result, int idx, int n)
{
	double s = 0; int i;
	int k;
	for (k = idx; k < n; k = k + 2)
	{
		for (i = 0; i < n; i++)
		{
			s = s + A[k][i] * S[i];
		}
		Result[k] = s;
		s = 0;
	}
}
int main()
{
	int i = 0; int j = 0; double k = 10; double eps = 0.001;
	double Norm1, Norm2;
	double **A, *B, *R, *X, *AZ, *Z, beta;
	double s3 = 0, s1 = 0, s2 = 0; int t = 1; 	double alpha;
	double *scalar;
	int n;
	cout << "Size" << endl;
	cin >> n;
	A = new double*[n];
	for (i = 0; i < n; i++)
		A[i] = new double[n];
	AZ = new double[n];
	Z = new double[n];
	X = new double[n];
	R = new double[n];
	B = new double[n];
	double *scalar1 = new double[1];
	double *Result = new double[n];
	for (i = 0; i < n; i++)
		Result[i] = 0;
	for (i = 0; i < n; i++)
		X[i] = 0;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i != j)
				A[i][j] = 1;
			else
				A[i][j] = 2 * (i + 1);
		}
		B[i] = t++;
	}
	

	int start_time = clock();
	thread t1(Func, A, X, Result, 0, n);
	thread t2(Func, A, X, Result, 1, n);

	t1.join();
	t2.join();
	for (i = 0; i < n; i++)
	{
		R[i] = B[i] - Result[i];
	}

	cout << endl << endl;
	for (i = 0; i < n; i++)
	{
		Z[i] = R[i];
	}
	scalar = new double[1];
	while (k >= eps)
	{

		scalar[0] = 0;
		thread p1(Func1, R, R, scalar, 0, n);
		thread p2(Func1, R, R, scalar, 1, n);

		p1.join();
		p2.join();

		thread t1(Func, A, Z, AZ, 0, n);
		thread t2(Func, A, Z, AZ, 1, n);
		t1.join();
		t2.join();

		scalar1[0] = 0;
		thread p11(Func1, AZ, Z, scalar1, 0, n);
		thread p21(Func1, AZ, Z, scalar1, 1, n);

		p11.join();
		p21.join();
		alpha = scalar[0] / scalar1[0];

		for (i = 0; i<n; i++)
		{
			X[i] = X[i] + alpha*Z[i];
			R[i] = R[i] - alpha*AZ[i];
		}

		for (i = 0; i < n; i++)
			Result[i] = 0;

		s1 = 0;

		scalar1[0] = 0;
		thread threadBeta1(Func1, R, R, scalar1, 0, n);
		thread threadBeta2(Func1, R, R, scalar1, 1, n);

		threadBeta1.join();
		threadBeta2.join();


		beta = scalar1[0] / scalar[0];

		for (i = 0; i<n; i++)

		{
			Z[i] = R[i] + beta*Z[i];
		}
		
		Norm1 = 0; Norm2 = 0;

		for (i = 0; i<n; i++)
		{
			Norm1 = Norm1 + R[i] * R[i];
			Norm2 = Norm2 + B[i] * B[i];
		}

		Norm1 = sqrt(Norm1);
		Norm2 = sqrt(Norm2);
		k = Norm1 / Norm2;

	}
	cout << "Answer is:" << endl;
	for (i = 0; i<n; i++)
	{
		cout << X[i] << " ";
		cout << endl;
	}

	cout << "time is: " << (clock() - start_time) / 1000.0 << endl;
	system("pause");
	return 0;
}
