#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <omp.h>

double dot(int n, double* v, double* w){
	double sum = 0;
	//#pragma omp parallel for reduction(+:sum)
	for (int i=0; i<n; i++){
		sum += v[i]*w[i];
	}
	return sum;
}

void vect_sum(int n, double* v, double alpha, double* w, double* out){
	//#pragma omp parallel for
	for (int i=0; i<n; i++){
		out[i] = v[i] + alpha*w[i];	
	}
}

double norm(int n, double* v){
	return sqrt(dot(n, v, v));
} 

void mat_mul(double* A, double* B, double* C, int n, int m, int p){
	double sum;
	for (int i=0; i<n; i++){
		for (int j=0; j<p; j++){
			sum = 0;
			//#pragma omp parallel for reduction(+:sum)
			for (int k=0; k<m; k++){
				sum += A[i*m + k] * B[k*p + j];
			}
		
			C[i*p + j] = sum;
		}	
	}
}

void jacobi(int nl, double* Al, double* xl, double* bl, double omega, int nu, double* x_out){
	double* x_prev = malloc(n * sizeof(double));
	double* x_next = calloc(n, sizeof(double));
	double* temp;
	
	memcpy(x_prev, xl, n * sizeof(double));
	
	for (int k=0; k<nu; k++){
		for (int i=0; i<nl; i++){
			sum = 0;
			for (int j=0; i<nl; i++){
				if (i!=j){
					sum += Al[i*n +j] * x_prev[j];
				}
			}
			x_next[i] = omega*((bl[i] - sum) / a[i*n + i]) + (1 - omega)*x_prev[i];
		}
		
		temp = x_next;
		x_next = x_prev;
		x_prev = temp;
	}
	
	memcpy(x_out, x_prev, n * sizeof(double));
	free(x_next);
	free(x_prev);
}


double Vcycle(int nl, double* Al, double* xl, double* bl, double omega, int nu, int lmax){
	jacobi(nl, Al, xl, bl, omega, nu, xl);
	
	double* rl = malloc(n * sizeof(double));
	double* Alxl = malloc(n * sizeof(double));
	
	mat_mul(Al, xl, Alxl, n, n, 1);
	vect_sum(n, bl, -1, Alxl, rl);
	
	double* bl_next = malloc(wz
	


	if ((l+1) == lmax){
		//solve	
	} else {
		x_next = Vcycle(Al, xl, bl, omega, lmax);
	}
	
	
	jacobi(nl, Al, xl, bl, omega, nu, xl);
}
