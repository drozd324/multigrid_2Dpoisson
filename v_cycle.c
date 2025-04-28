#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

double dot(int n, double* v, double* w){
	double sum = 0;
	for (int i=0; i<n; i++){
		sum += v[i]*w[i];
	}
	return sum;
}

void vect_sum(int n, double* v, double alpha, double* w, double* out){
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
 * @param int     n     Shape of square matrix Al
 * @param double* A     Matrix of size nl x nl
 * @param double* x     Initial guess solution to system of length nl
 * @param double* b     Vector of lenght nl
 * @param double  omega Parameter for weighted jacobi. For non weighted jacobi = 1
 * @param double  nu    Max number of iteraions
 * @param double  eps   Maximum resiual
 * @param double* x_out Output solution vector
*/
void jacobi(int n, double* A, double* x, double* b, double omega, int nu, double eps, double* x_out){
	double* x_prev = calloc(n, sizeof(double));
	double* x_next = calloc(n, sizeof(double));
	double* temp;
	double* residual = malloc(n * sizeof(double)); 
	
	memcpy(x_prev, x, n * sizeof(double));

	double sum; 
	for (int k=0; k<nu; k++){
		for (int i=0; i<n; i++){
			sum = 0;
			for (int j=0; i<n; i++){
				if (i!=j){
					sum += A[i*n +j] * x_prev[j];
				}
			}
			x_next[i] = omega*((b[i] - sum) / A[i*n + i]) + (1 - omega)*x_prev[i];
		}
		
		temp = x_next;
		x_next = x_prev;
		x_prev = temp;		
		
		// check residual
		mat_mul(A, x_prev, residual, n, n, 1);
		vect_sum(n, b, -1, residual, residual);
		if (norm(n, residual) < eps){
			break;
		}
		
	}
	
	memcpy(x_out, x_prev, n * sizeof(double));
	free(x_next);
	free(x_prev);
}


double f(double x1, double x2){
	return 2 * M_PI * M_PI * sin(M_PI * x1) * sin(M_PI * x2);  
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
			b[i*N + j] = f((double)(i/N), (double)(j/N));
		}
	}
}
//void restriction(int N, double* u, double* out){
//	int n = N/2; 
//	int i;
//	int j
//	int I;
//	int J;
//
//	// interior
//	for (i=1; i<(n-1); i++){
//		I = i*2;
//		for (j=1; j<(n-1); j++){
//			J = j*2;
//			out[i*n + j] = (1/18) * (u[N*(I+1) + (J-1)] +  u[N*(I+1) + J] + u[N*(I+1) + (J+1)] +
//									 u[N*(I  ) + (J-1)]+10*u[N*(I  ) + J] + u[N*(I  ) + (J+1)] +
//									 u[N*(I-1) + (J-1)] +  u[N*(I-1) + J] + u[N*(I-1) + (J+1)]);
//		}
//	}
//
//	// left
//	j = 0;
//	J = J*2;
//	for (i=1; i<(n-1); i++){
//		I = i*2;
//		out[i*n + j] = (1/12) * (						u[N*(I+1) + J] + u[N*(I+1) + (J+1)] +
//													  7*u[N*(I  ) + J] + u[N*(I  ) + (J+1)] +
//														u[N*(I-1) + J] + u[N*(I-1) + (J+1)]);
//	}
//
//	// right
//	j = n;
//	J = j*2;
//	for (i=1; i<(n-1); i++){
//		I = i*2;
//		out[i*n + j] = (1/12) * (u[N*(I+1) + (J-1)] + u[N*(I+1) + J] +
//							 	 u[N*(I  ) + (J-1)]+7*u[N*(I  ) + J] +
//								 u[N*(I-1) + (J-1)] + u[N*(I-1) + J] );
//	}
//
//	// top 
//	i = 0;
//	I = i*2;
//	for (j=1; j<(n-1); j++){
//		J = j*2;
//		out[i*n + j] = (1/12) * (
//								u[N*(I  ) + (J-1)] + u[N*(I  ) + J] + u[N*(I  ) + (J+1)] +
//								u[N*(I-1) + (J-1)]+7*u[N*(I-1) + J] + u[N*(I-1) + (J+1)]);
//	}
//
//	// bottom
//	i = n;
//	I = i*2;
//	for (j=1; j<(n-1); j++){
//		J = j*2;
//		out[i*n + j] = (1/12) * (u[N*(I+1) + (J-1)] + u[N*(I+1) + J] + u[N*(I+1) + (J+1)] +
//								u[N*(I  ) + (J-1)]+7*u[N*(I  ) + J] + u[N*(I  ) + (J+1)] 
//																						);
//	}
//
//
//	// top left
//	i = 0;
//	j = 0;
//	I = i*2;
//	J = j*2;
//	out[i*n + j] = (1/8) * (   
//										       5*u[N*(I  ) + J] + u[N*(I  ) + (J+1)] +
//											     u[N*(I-1) + J] + u[N*(I-1) + (J+1)]);
//
//	// bottom left
//	i = n;
//	j = 0;
//	I = i*2;
//	J = j*2;
//	out[i*n + j] = (1/8) * (				5*u[N*(I+1) + J] + u[N*(I+1) + (J+1)] +
//											  u[N*(I  ) + J] + u[N*(I  ) + (J+1)] +
//																				 );
//
//	// bottom right
//	i = n;
//	j = n;
//	I = i*2;
//	J = j*2;
//	out[i*n + j] = (1/8) * (u[N*(I+1) + (J-1)] + u[N*(I+1) + J] +
//							u[N*(I  ) + (J-1)]+5*u[N*(I  ) + J] 
//																					);
//
//	// top right
//	i = 0;
//	j = n;
//	I = i*2;
//	J = j*2;
//	out[i*n + j] = (1/8) * (
//							u[N*(I  ) + (J-1)]+5*u[N*(I  ) + J] +
//							u[N*(I-1) + (J-1)] + u[N*(I-1) + J]                 );
//
//}


void restriction(int N, double* u, double* out){
	int n = N/2; 
	int i;
	int j;
	int I;
	int J;

	for (i=0; i<n; i++){
		I = i*2;
		for (j=0; j<n; j++){
			J = j*2;
			out[i*n + j] = u[I*N + J];
		}
	}
}

void prolongate(int n, double* u, double* out){
	int N = n*2; 
	int i;
	int j;
	int I;
	int J;

	// instert
	for (i=0; i<N; i+=2){
		I = i/2;
		for (j=0; j<N; j+=2){
			J = j/2;
			out[i*N + j] = u[I*n + J]; 	
		}
	}

	// vertical
	for (i=0; i<(N-2); i+=2){
		I = i/2;
		for (j=0; j<N; j+=2){
			J = j/2;
			out[(i+1)*N + j] = (u[(I)*n + J] + u[(I+1)*n + J]) / 2;
		}
	}

	// horizontal
	for (i=0; i<N; i+=2){
		I = i/2;
		for (j=0; j<(N-2); j+=2){
			J = j/2;
			out[(i+1)*N + j] = (u[I*n + J] + u[I*n + (J+1)]) / 2;
		}
	}

	// middle
	for (i=0; i<(N-2); i+=2){
		I = i/2;
		for (j=0; j<(N-2); j+=2){
			J = j/2;
			out[(i+1)*N + (j+1)] = (u[(I  )*n + (J  )] + u[(I  )*n + (J+1)] +
									u[(I+1)*n + (J  )] + u[(I+1)*n + (J+1)]) / 4; 	
		}
	}
}


void Vcycle(int nl, double* Al, double* xl, double* bl, 
			double omega, int nu, int lmax, int l, double eps){
	int nl_next = nl/2;
	int size = nl*nl;
	int size_next = nl_next*nl_next;

	double* Al_next = malloc(size_next*size_next * sizeof(double));
	double* bl_next = malloc(size_next* sizeof(double));
	double* xl_next = malloc(size_next* sizeof(double));
	double* rl = malloc(size * sizeof(double));

	make_matrix(Al_next, bl_next, nl_next);
	
	jacobi(size, Al, xl, bl, omega, nu, -1, xl);
	
	mat_mul(Al, xl, rl, size, size, 1);
	vect_sum(size, bl, -1, rl, rl);
	
	restriction(size, rl, bl_next);                        
	
	if ((l+1) == lmax){
		jacobi(size_next, Al_next, xl_next, bl_next, omega, nu, eps, xl_next);
	} else {
		Vcycle(size_next, Al_next, xl_next, bl_next, omega, nu, lmax, l+1, eps);
	}
	
	double* Pxl_next = malloc(size_next* sizeof(double));
	prolongate(size, xl_next, Pxl_next);
	vect_sum(size, xl, 1, Pxl_next, xl);
	
	jacobi(size, Al, xl, bl, omega, nu, -1, xl);


	free(rl);
	free(bl_next);
	free(xl_next);
	free(Al_next);
	free(Pxl_next);
}
