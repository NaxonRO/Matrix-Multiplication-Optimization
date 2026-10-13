/*
 * Tema 2 ASC
 * 2026 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
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
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				C[i * N + j] += A[k * N + i] * B[k * N + j];
			}
		}
	}

	// D = C * C^T
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			double sum = 0.0;
			for (int k = 0; k < N; k++) {
				sum += C[i * N + k] * C[j * N + k];
			}
			D[i * N + j] = sum;
			D[j * N + i] = sum;
		}
	}

	// y += D * row_i(C)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				y[j] += D[j * N + k] * C[i * N + k];
			}
		}
	}

	// y = y + x
	for (int i = 0; i < N; i++) {
		y[i] += x[i];
	}

	free(C);
	free(D);

	return y;
}
