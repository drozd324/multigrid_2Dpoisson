#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "v_cycle.h"

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
	int N = 128;
	
	FILE *fp = fopen("writeup/q2out.csv", "w");
	fprintf(fp, "N,time,num iter\n");
	
	for (int l=0; l<lmax; k++){

		long long int size = N*N;
		double* A = calloc(size * size, sizeof(double));
		double* u = calloc(size, sizeof(double));
		double* b = calloc(size, sizeof(double));
	
		make_matrix(A, b, N);
		
		double t1 = walltime();
		double Vcycle(N, A, x, b, omega, nu, lmax, EPS);
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
}

