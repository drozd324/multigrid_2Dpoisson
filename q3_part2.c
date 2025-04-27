#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "conj_grad.h"

double f(double x1, double x2){
	return 2 * M_PI * M_PI * sin(M_PI * x1) * sin(M_PI * x2);  
}
	
void print_mat(int n, double* A){
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			printf("%lf, ", A[i*n + j]);
		}
		printf("\n");
	}
}

#define GRID_PTS 6
#define EPS 1e-7

int main(){
	int N = 2;
	
	FILE *fp = fopen("writeup/q2out.csv", "w");
	FILE *fp_sol = fopen("writeup/q2sol.txt", "w");
	fprintf(fp, "N,time,num iter\n");
	
	for (int k=0; k<GRID_PTS; k++){
		N *= 2;

		long long int size = N*N;
		double* A = calloc(size * size, sizeof(double));
		double* u = calloc(size, sizeof(double));
		double* b = calloc(size, sizeof(double));

		for (int i=0; i<size; ++i){
			A[i*size + i] =  4.0;
		}
		for (int i=0; i<size-1; ++i){
			if ((i % N) != N-1){
				A[(i+1)*size + i+0] = 1.0;
				A[(i+0)*size + i+1] = 1.0;	
			}
		}
		for (int i=0; i<size-N; ++i){
			A[(i+N)*size + i+0] = 1.0;
			A[(i+0)*size    + i+N] = 1.0;
		}

		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				b[i*N + j] = f((double)i/N, (double)j/N);
			}
		}

		
		int num_iter;
		double t1 = walltime();
		conjugate_gradient(size, A, b, u, EPS, u, &num_iter);
		double time = walltime() - t1;
	
		printf("N = %d Converged at iter = %d\n", N, num_iter);

		if (k == GRID_PTS-1){	
			for (int i=0; i<size; i++){
				fprintf(fp_sol, "%lf ", u[i]);
			}	
			fprintf(fp_sol, "\n");
		}

		free(A);
		free(u);
		free(b);
		fprintf(fp, "%d,%lf,%d\n", N, time, num_iter);
	}
	fclose(fp);
	fclose(fp_sol);
}

