/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>
#include <stdlib.h>

/*
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B, double *x) {
	double *C = calloc(N * N, sizeof(double));
	double *D = calloc(N * N, sizeof(double));
	double *y = calloc(N, sizeof(double));

	if (!C || !D || !y) {
		free(C);
		free(D);
		free(y);
		return NULL;
	}

	// C = A^T * B
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, A, N, B, N, 0.0, C, N);

	// D = C * C^T
	cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, N, N, 1.0, C, N, 0.0, D, N);

	// sum_rows = sum_i row_i(C) = C^T * ones
	double *ones = malloc(N * sizeof(double));
	double *sum_rows = calloc(N, sizeof(double));
	if (!ones || !sum_rows) {
		free(ones);
		free(sum_rows);
		free(C);
		free(D);
		free(y);
		return NULL;
	}
	for (int i = 0; i < N; i++) ones[i] = 1.0;
	cblas_dgemv(CblasRowMajor, CblasTrans, N, N, 1.0, C, N, ones, 1, 0.0, sum_rows, 1);

	// y += D * sum_rows
	cblas_dsymv(CblasRowMajor, CblasUpper, N, 1.0, D, N, sum_rows, 1, 1.0, y, 1);

	// y = y + x
	cblas_daxpy(N, 1.0, x, 1, y, 1);

	free(ones);
	free(sum_rows);
	free(C);
	free(D);

	return y;
}
