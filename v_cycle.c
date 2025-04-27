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

/*
 * @brief Weighted Jacobi iteration function
 *
 * @param int     nl    Shape of square matrix Al
 * @param double* Al    Matrix of size nl x nl
 * @param double* xl    Initial guess solution to system of length nl
 * @param double* bl    Vector of lenght nl
 * @param double  omega Parameter for weighted jacobi. For non weighted jacobi = 1
 * @param double  nu    Max number of iteraions
 * @param double  r     Maximum resiual
 * @param double* x_out Output solution vector
*/
void jacobi(int nl, double* Al, double* xl, double* bl, double omega, int nu, double eps, double* x_out){
	double* x_prev = calloc(nl, sizeof(double));
	double* x_next = calloc(nl, sizeof(double));
	double* temp;
	double* residual = malloc(nl * sizeof(double)); 
	
	memcpy(x_prev, xl, n * sizeof(double));
	
	for (int k=0; k<nu; k++){
		for (int i=0; i<nl; i++){
			sum = 0;
			for (int j=0; i<nl; i++){
				if (i!=j){
					sum += Al[i*n +j] * x_prev[j];
				}
			}
			x_next[i] = omega*((bl[i] - sum) / Al[i*n + i]) + (1 - omega)*x_prev[i];
		}
		
		temp = x_next;
		x_next = x_prev;
		x_prev = temp;		
		
		// check residual
		mat_mul(Al, x_prev, residual, nl, nl, 1);
		vect_sum(nl, bl, -1, residual, residual);
		if (norm(resiual) < eps){
			break;
		}
		
	}
	
	memcpy(x_out, x_prev, n * sizeof(double));
	free(x_next);
	free(x_prev);
}

/*
 * @brief Function to construct matrix for 2d poisson problem
 *
 * @param int     N    Shape of square matrix Al
 * @param double* A    Matrix of size N x N 
 * @param double* b    Vector of lenght N
*/
void make_matrix(double* A, double* b, int N){
	int size = N*N;

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
}

void restriction(int N, double* u, double* out){
	int n = N/2; 
	int I;
	int J;

	for (int i=1; i<(n-1); i++){
		I = i*2;
		for (int j=1; j<(n-1); j++){
			J = j*2;
			out[i*n + j] = (1/9) * (u[N*(I+1) + (J-1)] + u[N*(I+1) + J] + u[N*(I+1) + (J+1)] +
									u[N*(I  ) + (J-1)] + u[N*(I  ) + J] + u[N*(I  ) + (J+1)] +
									u[N*(I-1) + (J-1)] + u[N*(I-1) + J] + u[N*(I-1) + (J+1)]);
		}
	}

	i = 0
	I = i*2;
	for (int j=1; j<(n-1); j++){
		J = j*2;
		out[i*n + j] = (1/6) * (						u[N*(I+1) + J] + u[N*(I+1) + (J+1)] +
														u[N*(I  ) + J] + u[N*(I  ) + (J+1)] +
														u[N*(I-1) + J] + u[N*(I-1) + (J+1)]);
	}

	i = n
	I = i*2;
	for (int j=1; j<(n-1); j++){
		J = j*2;
		out[i*n + j] = (1/6) * (u[N*(I+1) + (J-1)] + u[N*(I+1) + J] +
								u[N*(I  ) + (J-1)] + u[N*(I  ) + J] +
								u[N*(I-1) + (J-1)] + u[N*(I-1) + J] );
	}

}

void Vcycle(int nl, double* Al, double* xl, double* bl, double omega, int nu, int lmax, double r, double eps){
	jacobi(nl, Al, xl, bl, omega, nu, eps, xl);
	
	double* rl = malloc(n * sizeof(double));
	
	mat_mul(Al, xl, rl, n, n, 1);
	vect_sum(n, bl, -1, rl, rl);
	
	double* bl_next = malloc(((nl/2) + 1) * sizeof(double));
	
	if ((l+1) == lmax){
		//solve	
	} else {
		x_next = Vcycle(Al, xl, bl, omega, lmax);
	}
	
	
	jacobi(nl, Al, xl, bl, omega, nu, xl);
}
