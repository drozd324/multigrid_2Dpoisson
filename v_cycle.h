#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <omp.h>

#define MAX_ITER 10000

double dot(int n, double* v, double* w);
void vect_sum(int n, double* v, double alpha, double* w, double* out);
double norm(int n, double* v);
void mat_mul(double* A, double* B, double* C, int n, int m, int p);
int jacobi(int nl, double* Al, double* xl, double* bl, double omega, int nu, double eps, double* x_out);
void make_matrix(double* A, double* b, int N);
int Vcycle(int nl, double* Al, double* xl, double* bl, 
			double omega, int nu, int lmax, int l, double eps);
