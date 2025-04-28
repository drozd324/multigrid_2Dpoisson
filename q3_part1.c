#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "v_cycle.h"

	
void print_mat(int n, double* A){
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			printf("%lf, ", A[i*n + j]);
		}
		printf("\n");
	}
}


/**
 * @brief Function to measure time
 *
 * @param[out] Actual time in seconds
 */
double walltime(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double wtime = (double) (t.tv_sec + t.tv_usec*1e-6);
    return wtime;
}

#define GRID_PTS 6
#define EPS 1e-7

int main(){
	int N = 64;
	int lmax = 5;
	double omega = 1;
	double nu = 4;
	
	//FILE *fp = fopen("writeup/q2out.csv", "w");
	//fprintf(fp, "N,time,num iter\n");
	
	for (int l=2; l<lmax; l++){

		long long int size = N*N;
		double* A = calloc(size * size, sizeof(double));
		double* u = calloc(size, sizeof(double));
		double* b = calloc(size, sizeof(double));
		
		double t1 = walltime();
		make_matrix(A, b, N);
		Vcycle(N, A, u, b, omega, nu, lmax, l, EPS);
		//Vcycle(int nl, double* Al, double* xl, double* bl, 
		//	double omega, int nu, int lmax, int l, double r, double eps);
		double time = walltime() - t1;
	
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

