/*
 * Tema 2 ASC
 * 2026 Spring
 */
#include "utils.h"
#include <stdlib.h>

/*
 * Add your optimized implementation here
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
	double *Ak = A, *Bk = B;
	for (int k = 0; k < N; k++, Ak += N, Bk += N) {
		double *Ci = C;
		for (int i = 0; i < N; i++, Ci += N) {
			double a = Ak[i];
			double *pci = Ci, *pbk = Bk;
			int j = 0;
			for (; j + 7 < N; j += 8, pci += 8, pbk += 8) {
				pci[0] += a * pbk[0]; pci[1] += a * pbk[1];
				pci[2] += a * pbk[2]; pci[3] += a * pbk[3];
				pci[4] += a * pbk[4]; pci[5] += a * pbk[5];
				pci[6] += a * pbk[6]; pci[7] += a * pbk[7];
			}
			for (; j < N; j++, pci++, pbk++) {
				*pci += a * *pbk;
			}
		}
	}

	// D = C * C^T
	for (int i = 0; i < N; i++) {
		double *Ci = C + i * N;
		int j = i;
		for (; j + 7 < N; j += 8) {
			double s0=0, s1=0, s2=0, s3=0, s4=0, s5=0, s6=0, s7=0;
			double *pci = Ci;
			double *pj0 = C + j*N,     *pj1 = C + (j+1)*N;
			double *pj2 = C + (j+2)*N, *pj3 = C + (j+3)*N;
			double *pj4 = C + (j+4)*N, *pj5 = C + (j+5)*N;
			double *pj6 = C + (j+6)*N, *pj7 = C + (j+7)*N;
			int k = 0;
			for (; k + 3 < N; k += 4,
				pci += 4, pj0 += 4, pj1 += 4, pj2 += 4, pj3 += 4,
				pj4 += 4, pj5 += 4, pj6 += 4, pj7 += 4) {
				double c0=pci[0], c1=pci[1], c2=pci[2], c3=pci[3];
				s0 += c0*pj0[0] + c1*pj0[1] + c2*pj0[2] + c3*pj0[3];
				s1 += c0*pj1[0] + c1*pj1[1] + c2*pj1[2] + c3*pj1[3];
				s2 += c0*pj2[0] + c1*pj2[1] + c2*pj2[2] + c3*pj2[3];
				s3 += c0*pj3[0] + c1*pj3[1] + c2*pj3[2] + c3*pj3[3];
				s4 += c0*pj4[0] + c1*pj4[1] + c2*pj4[2] + c3*pj4[3];
				s5 += c0*pj5[0] + c1*pj5[1] + c2*pj5[2] + c3*pj5[3];
				s6 += c0*pj6[0] + c1*pj6[1] + c2*pj6[2] + c3*pj6[3];
				s7 += c0*pj7[0] + c1*pj7[1] + c2*pj7[2] + c3*pj7[3];
			}
			for (; k < N; k++, pci++, pj0++, pj1++, pj2++, pj3++,
				pj4++, pj5++, pj6++, pj7++) {
				double cik = *pci;
				s0 += cik * *pj0; s1 += cik * *pj1;
				s2 += cik * *pj2; s3 += cik * *pj3;
				s4 += cik * *pj4; s5 += cik * *pj5;
				s6 += cik * *pj6; s7 += cik * *pj7;
			}
			D[i*N+j]   = s0; D[j*N+i]     = s0;
			D[i*N+j+1] = s1; D[(j+1)*N+i] = s1;
			D[i*N+j+2] = s2; D[(j+2)*N+i] = s2;
			D[i*N+j+3] = s3; D[(j+3)*N+i] = s3;
			D[i*N+j+4] = s4; D[(j+4)*N+i] = s4;
			D[i*N+j+5] = s5; D[(j+5)*N+i] = s5;
			D[i*N+j+6] = s6; D[(j+6)*N+i] = s6;
			D[i*N+j+7] = s7; D[(j+7)*N+i] = s7;
		}
		for (; j + 3 < N; j += 4) {
			double s0=0, s1=0, s2=0, s3=0;
			double *pci = Ci;
			double *pj0 = C + j*N,     *pj1 = C + (j+1)*N;
			double *pj2 = C + (j+2)*N, *pj3 = C + (j+3)*N;
			int k = 0;
			for (; k + 3 < N; k += 4,
				pci += 4, pj0 += 4, pj1 += 4, pj2 += 4, pj3 += 4) {
				double c0=pci[0], c1=pci[1], c2=pci[2], c3=pci[3];
				s0 += c0*pj0[0] + c1*pj0[1] + c2*pj0[2] + c3*pj0[3];
				s1 += c0*pj1[0] + c1*pj1[1] + c2*pj1[2] + c3*pj1[3];
				s2 += c0*pj2[0] + c1*pj2[1] + c2*pj2[2] + c3*pj2[3];
				s3 += c0*pj3[0] + c1*pj3[1] + c2*pj3[2] + c3*pj3[3];
			}
			for (; k < N; k++, pci++, pj0++, pj1++, pj2++, pj3++) {
				double cik = *pci;
				s0 += cik * *pj0; s1 += cik * *pj1;
				s2 += cik * *pj2; s3 += cik * *pj3;
			}
			D[i*N+j]   = s0; D[j*N+i]     = s0;
			D[i*N+j+1] = s1; D[(j+1)*N+i] = s1;
			D[i*N+j+2] = s2; D[(j+2)*N+i] = s2;
			D[i*N+j+3] = s3; D[(j+3)*N+i] = s3;
		}
		for (; j < N; j++) {
			double sum = 0.0;
			double *pci = Ci, *pcj = C + j * N;
			int k = 0;
			for (; k + 3 < N; k += 4, pci += 4, pcj += 4)
				sum += pci[0]*pcj[0] + pci[1]*pcj[1] + pci[2]*pcj[2] + pci[3]*pcj[3];
			for (; k < N; k++, pci++, pcj++)
				sum += *pci * *pcj;
			D[i*N+j] = sum;
			D[j*N+i] = sum;
		}
	}

	// y += D_sym * row_i(C)
	for (int i = 0; i < N; i++) {
		double *v = C + i * N;
		for (int j = 0; j < N; j++) {
			double *dj  = D + j * N;
			double  vj  = v[j];
			double  tmp = dj[j] * vj;
			int k = j + 1;
			for (; k + 7 < N; k += 8) {
				double d0=dj[k],   d1=dj[k+1], d2=dj[k+2], d3=dj[k+3];
				double d4=dj[k+4], d5=dj[k+5], d6=dj[k+6], d7=dj[k+7];
				double v0=v[k],    v1=v[k+1],  v2=v[k+2],  v3=v[k+3];
				double v4=v[k+4],  v5=v[k+5],  v6=v[k+6],  v7=v[k+7];
				tmp    += d0*v0 + d1*v1 + d2*v2 + d3*v3
				        + d4*v4 + d5*v5 + d6*v6 + d7*v7;
				y[k]   += d0 * vj;
				y[k+1] += d1 * vj;
				y[k+2] += d2 * vj;
				y[k+3] += d3 * vj;
				y[k+4] += d4 * vj;
				y[k+5] += d5 * vj;
				y[k+6] += d6 * vj;
				y[k+7] += d7 * vj;
			}
			for (; k < N; k++) {
				double dval = dj[k];
				tmp  += dval * v[k];
				y[k] += dval * vj;
			}
			y[j] += tmp;
		}
	}

	// y = y + x
	for (int i = 0; i < N; i++)
		y[i] += x[i];

	free(C);
	free(D);

	return y;
}
