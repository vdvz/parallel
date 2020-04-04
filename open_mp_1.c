#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
 
#define N 4000


int main(int argc, char* argv[]) {
	int i,j;
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

	int count = 0;
	double sum_B = 0.0;

	struct timeval start, end;
	gettimeofday(&start, NULL);

#pragma omp parallel private(i,j) 
	

#pragma omp for 
	for (i = 0; i < N; i++)
		sum_B += B[i] * B[i];
	double sub_B = sqrt(sum_B);
	

#pragma omp for
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			x[i] += A[i][j] * X[j];
		}
	}
	

#pragma omp for
	for (i = 0; i < N; i++)
		x[i] = x[i] - B[i];
	
	double sum_x = 0;
#pragma omp for
	for (i = 0; i < N; i++)
		sum_x += x[i] * x[i];
	double sub_x = sqrt(sum_x);
	double c_div = sub_x / sub_B;
	double p_div = c_div;
	while (c_div >= eps) {
		count++;
		//mult_matrix_btvector
#pragma omp for
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				x[i] += A[i][j] * X[j];
			}
		}
		//sub_vector
#pragma omp for
		for (i = 0; i < N; i++)
			x[i] = x[i] - B[i];
		//mult_vector_bynumber
#pragma omp for
		for (i = 0; i < N; i++)
			x[i] *= t;
		//sub__vector
#pragma omp for
		for (i = 0; i < N; i++)
			X[i] = X[i] - x[i];
		//find_solution
		if (c_div > p_div) t = -0.01;
		p_div = c_div;
		//check 
		//mult_matrix_byvector
#pragma omp for
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				x[i] += A[i][j] * X[j];
			}
		}
		//sub_vector(x, B, x);
#pragma omp for
		for (i = 0; i < N; i++)
			x[i] = x[i] - B[i];
		//mod x
		double sum_x = 0;
#pragma omp for
		for (i = 0; i < N; i++)
			sum_x += x[i] * x[i];
		double sub_x = sqrt(sum_x);
		c_div = sub_x / sub_B;
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