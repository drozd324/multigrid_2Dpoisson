#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h> 
#include <unistd.h>
#include "v_cycle.h"

#define EPS 1e-7

int main(int argc, char *argv[]) {

	int N = 16;
	int lmax = 2;
	double omega = 1;
	double nu = 4;
	int option;
  
    while ((option = getopt(argc, argv, "L:")) != -1) {
        switch (option) {
            case 'L':
	            lmax = atoi(optarg);
				break;
		}
	}
	
	for (int i=0; i<4; i++){
		
		printf("\nN = %d\n", N);
	
		long long int size = N*N;
		double* A = calloc(size * size, sizeof(double));
		double* u = calloc(size, sizeof(double));
		double* b = calloc(size, sizeof(double));
		make_matrix(A, b, N);

		double t0 = clock();
		Vcycle(N, A, u, b, omega, nu, lmax, 0, EPS);
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

		N *= 2;
	}
}

