#include "/home/andrew/OpenBLAS/include/cblas.h"

void square_dgemm(int n, double *A, double *B, double *C)
{
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, A, n, B, n, 1, C, n);
}
