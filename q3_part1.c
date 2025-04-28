#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h> 
#include "v_cycle.h"

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
	int N = 64;//128; // 2^7
	int lmax = 6;
	double omega = 1;
	double nu = 4;
	
	//FILE *fp = fopen("writeup/q2out.csv", "w");
	//fprintf(fp, "N,time,num iter\n");
	
	for (int l=2; l<lmax; l++){

		long long int size = N*N;
		double* A = calloc(size * size, sizeof(double));
		double* u = calloc(size, sizeof(double));
		double* b = calloc(size, sizeof(double));
		make_matrix(A, b, N);

		double t0 = clock();
		Vcycle(N, A, u, b, omega, nu, l, 0, EPS);
		double t1 = clock();  
		double time = (t1 - t0) / (CLOCKS_PER_SEC);
	
		printf("%lf\n", time);

//		for (int i=0; i<size; i++){
//			fprintf(fp_sol, "%lf ", u[i]);
//		}	
//		fprintf(fp_sol, "\n");

		free(A);
		free(u);
		free(b);
		//fprintf(fp, "%d,%lf,%d\n", N, time, num_iter);
	}
	//fclose(fp);
}

