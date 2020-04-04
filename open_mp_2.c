#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define N 4000

void mult_matrix_vector(double** A, double* B, double* vector);

void mult_vector_number(double* vector, double num);

void sub_vector(double* vector_A, double* B, double* vector);

void result(double** A, double* X, double* B, double* x, double t);

double check(double** A, double* X, double* B, double mod_vector_B);

double mod(double* vector);

void sub_vector(double* vector_A, double* B, double* vector);

int main(int argc, char* argv[]) {
	int i, j;
	double t = 0.00001;
	double eps = 0.00001;
	
	double A[N][N];
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i == j) A[i][j] = 2.0;
			else A[i][j] = 1.0;
		}
	}

	double X[N];
	for (i = 0; i < N; i++) {
		X[i] = 0.0;
	}

	double B[N];
	for (i = 0; i < N; i++) {
		B[i] = N + 1.0;
	}

	double x[N];
	for (i = 0; i < N; i++) {
		x[i] = 0.0;
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);

	double mod_vector_B = mod(B);
	double c_div = check(A, X, B, mod_vector_B);
	double p_div = c_div;
	int count = 0;
	while (c_div >= eps) {
		count++;
		result(A, X, B, x, t);
		if (c_div > p_div) t = -0.01;
		p_div = c_div;
		c_div = check(A, X, B, mod_vector_B);
	}

	gettimeofday(&end, NULL);
	double dt_sec = (end.tv_sec - start.tv_sec);
	double dt_usec = (end.tv_usec - start.tv_usec);
	double dt = dt_sec + 0.000001 * dt_usec;
	printf("time diff %lf \n", dt);

	printf("%d\n", count);
	for (i = 0; i < N; i++)
		printf("%lf ", X[i]);

	return 0;

}

void mult_matrix_vector(double** A, double* B, double* vector) {
	int i, j;
#pragma omp parallel for schedule(static, (int)N/omp_get_num_threads()) private(j)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			vector[i] += A[i][j] * B[j];
		}
	}
}
void mult_vector_number(double* vector, double num) {
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < N; i++)
		vector[i] *= num;
}
void sub_vector(double* vector_A, double* B, double* vector) {
	int i;
#pragma omp parallel for private(i) 
	for (i = 0; i < N; i++)
		vector[i] = vector_A[i] - B[i];
}
double mod(double* vector) {
	double sum = 0;
	int i;
#pragma omp parallel for reduction(+:sum)
	for (i = 0; i < N; i++)
		sum += vector[i] * vector[i];
	return sqrt(sum);
}
double check(double** A, double* X, double* B, double mod_vector_B) {
	double* x = (double*)calloc(N, sizeof(double));
	mult_matrix_vector(A, X, x);
	sub_vector(x, B, x);
	double mod_tmp = mod(x);
	free(x);
	return mod_tmp / mod_vector_B;
}

void result(double** A, double* X, double* B, double* x, double t) {
	mult_matrix_vector(A, X, x);
	sub_vector(x, B, x);
	mult_vector_number(x, t);
	sub_vector(X, x, X);
}
