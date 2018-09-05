
/*----------------------------------------------------------------*/

#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <omp.h>

/*-----------------------------------------------------------*/

#include "hclfftw.h"
#include "initmatrix.h"
#include "transpose.h"
#include "fpm.h"

/*-----------------------------------------------------------*/

int main(int argc, char** argv)
{
    if (argc != 6)
    {
       fprintf(
          stderr,
          "Usage: %s <n> <ngroups> <nthreadspergroup> "
          "<verbosity (0|1)> "
          "Number of invocations.\n",
          argv[0]);
       exit(EXIT_FAILURE);
    }

    int N = atoi(argv[1]);
    unsigned int nGroups = atoi(argv[2]);
    unsigned int nThreadsPerGroup = atoi(argv[3]);
    unsigned int verbosity = atoi(argv[4]);
    unsigned int numInvocations = atoi(argv[5]);

    fftw_complex *lMatrix = (fftw_complex*)fftw_malloc(
                           N*N*sizeof(fftw_complex));

    if (lMatrix == NULL)
    {
       fprintf(
         stderr,
         "Memory allocation failure.\n"
       );
       exit(EXIT_FAILURE);
    }

    fftw_init_threads();
    fftw_plan_with_nthreads(nThreadsPerGroup);

    hclFillSignal2D(
       N, N, nThreadsPerGroup*nGroups, lMatrix);

    int* rowd = (int*)malloc(nGroups*sizeof(int));
    determineRowDistributionLB(
       N, nGroups, rowd);

    double* tGroups = (double*)calloc(nGroups, sizeof(double));

    double tstart, tend;
    struct timeval start, end;

    gettimeofday(&start, NULL);

    int invocation;
    for (invocation  = 0; invocation < numInvocations; invocation++)
    {
       fftwlocal(
          FFTW_FORWARD, rowd, N,
          nThreadsPerGroup, nGroups, lMatrix, tGroups);
    }

    gettimeofday(&end, NULL);

    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;

    fftw_cleanup_threads();
    free(rowd);
    fftw_free(lMatrix);
    free(tGroups);

    double t = (tend - tstart) / (double)numInvocations;
    double pSize = (double)N*(double)N;
    double speed = 2.5 * pSize * log2(pSize) / t;
    double mflops = speed * 1e-06;

    printf(
        "n=%d, time(sec)=%0.6f, speed(mflops)=%0.6f\n",
        N, t, mflops
    );

    exit(EXIT_SUCCESS);
}

/*----------------------------------------------------------------*/
