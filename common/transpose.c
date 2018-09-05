
/*----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

/*----------------------------------------------------------------*/

#include "transpose.h"

/*----------------------------------------------------------------*/

void hcl_local_transpose_scalar_block(
    fftw_complex* X1,
    fftw_complex* X2,
    const int i, const int j,
    const int n,
    const int block_size,
    const unsigned int verbosity)
{
    int p, q;

    for (p = 0; p < min(n-i,block_size); p++) {
        for (q = 0; q < min(n-j,block_size); q++) {

           if (verbosity)
              printf(
                  "%d: i %d, j %d, p %d, q %d, index1 %d index2 %d\n",
                  omp_get_thread_num(),
                  i, j,
                  p, q,
                  i*n+j + p*n+q, j*n+i + q*n+p);

           double tmpr = X1[p*n+q][0];
           double tmpi = X1[p*n+q][1];
           X1[p*n+q][0] = X2[q*n+p][0];
           X1[p*n+q][1] = X2[q*n+p][1];
           X2[q*n+p][0] = tmpr;
           X2[q*n+p][1] = tmpi;
        }
    }
}

void hcl_transpose_block(
    fftw_complex* X,
    const int start, const int end,
    const int n,
    const unsigned int nt,
    const int block_size,
    const unsigned int verbosity)
{
    int i, j;

#pragma omp parallel for shared(X) private(i, j) num_threads(nt)
    for (i = 0; i < end; i += block_size) {
        for (j = 0; j < end; j += block_size) {
            if (verbosity)
               printf(
                  "%d: i %d, j %d\n",
                  omp_get_thread_num(), i, j);

            hcl_local_transpose_scalar_block(
                &X[start + i*n + j],
                &X[start + j*n + i], i, j, n, block_size, verbosity);
        }
    }
}

/*----------------------------------------------------------------*/
