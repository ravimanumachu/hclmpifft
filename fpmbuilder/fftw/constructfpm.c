
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
          "Usage: %s <n1> <n2> <ngroups> <nthreadspergroup>\n"
          "<verbosity (0|1)>\n"
          "Number of invocations.\n",
          argv[0]);
       exit(EXIT_FAILURE);
    }

    int N1 = atoi(argv[1]);
    int N2 = atoi(argv[2]);
    unsigned int nGroups = atoi(argv[3]);
    unsigned int nThreadsPerGroup = atoi(argv[4]);
    unsigned int verbosity = atoi(argv[5]);
    unsigned int numInvocations = atoi(argv[6]);

    fftw_complex *lMatrix = (fftw_complex*)fftw_malloc(
                           N1*N2*sizeof(fftw_complex));

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
       N1, N2, nThreadsPerGroup*nGroups, lMatrix);

    double* tGroups = (double*)calloc(nGroups, sizeof(double));
    int group;
    int* rowd = (int*)malloc(sizeof(int)*nGroups);
    for (group = 0; group < nGroups; group++)
    {
        rowd[group] = N1;
    }

    int invocation;
    for (invocation  = 0; invocation < numInvocations; invocation++)
    {
       fftwlocal(
          FFTW_FORWARD, rowd, N2,
          nThreadsPerGroup, nGroups, lMatrix, tGroups);
    }

    fftw_cleanup_threads();
    free(rowd);
    fftw_free(lMatrix);
    free(tGroups);

    unsigned int team;
    double avg = 0.;
    for (team = 0; team < nGroups; team++)
    {
        double t = tGroups[team] / (double)numInvocations;
        avg += t;

        double pSize = (double)N1*(double)N2;
        double speed = 2.5 * pSize * log2(pSize) / t;
        double mflops = speed * 1e-06;
        printf(
            "m=%d, n=%d, time(sec)=%0.6f, Speed(mflops)=%0.6f\n",
            N1, N2, t, mflops
        );
    }

    avg = avg / 4.0;
    double pSize = (double)N1*(double)N2;
    double speed = 2.5 * pSize * log2(pSize) / avg;
    double mflops = speed * 1e-06;

    printf(
        "m=%d, n=%d, time(sec)=%0.6f, Average speed(mflops)=%0.6f\n",
        N1, N2, avg, mflops
    );

    exit(EXIT_SUCCESS);
}

/*----------------------------------------------------------------*/
