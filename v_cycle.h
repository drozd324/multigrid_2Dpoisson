#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <omp.h>

double dot(int n, double* v, double* w);
void vect_sum(int n, double* v, double alpha, double* w, double* out);
double norm(int n, double* v);
void mat_mul(double* A, double* B, double* C, int n, int m, int p);
void jacobi(int nl, double* Al, double* xl, double* bl, double omega, int nu, double eps, double* x_out);
void make_matrix(double* A, double* b, int N);
double Vcycle(int nl, double* Al, double* xl, double* bl, double omega, int nu, int lmax);
