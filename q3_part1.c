#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h> 
#include "v_cycle.h"

#define EPS 1e-7

int main(){
	int N = 128;//128; // 2^7
	int lmax = 7;
	double omega = 1;
	double nu = 4;
	
	for (int l=2; l<lmax; l++){
	
		printf("\nl = %d\n", l);
	
		int size = N*N;
		double* A = calloc(size * size, sizeof(double));
		double* u = calloc(size, sizeof(double));
		double* b = calloc(size, sizeof(double));
		make_matrix(A, b, N);

		double t0 = clock();
		Vcycle(N, A, u, b, omega, nu, l, 0, EPS);
		double t1 = clock();
		double time = (t1 - t0) / (CLOCKS_PER_SEC);
	
		printf("time taken = %lf\n", time);
		
		double* r = calloc(size, sizeof(double));
		mat_mul(A, u, r, size, size, 1);
		vect_sum(size, b, -1, r, r);
		printf("final residual ||r|| = %lf", norm(size, r));
		printf(", ||r||/ len(r) = %lf\n", norm(size, r)/size);

		free(A);
		free(u);
		free(b);
		free(r);
	}
}

