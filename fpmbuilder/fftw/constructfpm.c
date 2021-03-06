
/*----------------------------------------------------------------*/

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
    if (argc != 7)
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
    unsigned int* rowd = (unsigned int*)malloc(sizeof(unsigned int)*nGroups);
    for (group = 0; group < nGroups; group++)
    {
        rowd[group] = N1 / nGroups;
    }

    double tstart, tend;
    struct timeval start, end;

    gettimeofday(&start, NULL);

    int invocation;
    for (invocation  = 0; invocation < numInvocations; invocation++)
    {
       fftwlocal(
          FFTW_FORWARD, rowd, N2,
          nThreadsPerGroup, nGroups, lMatrix, tGroups);
    }

    gettimeofday(&end, NULL);

    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;

    fftw_cleanup_threads();
    fftw_free(lMatrix);

    unsigned int team;
    double max = 0.;
    for (team = 0; team < nGroups; team++)
    {
        double t = tGroups[team] / (double)numInvocations;

        if (t > max)
        {
           max = t;
        }

        double pSize = (double)rowd[team]*(double)N2;
        double speed = 5.0 * pSize * log2(pSize) / t;
        double mflops = speed * 1e-06;
        printf(
            "Team %d: m=%d, n=%d, time(sec)=%0.6f, Speed(mflops)=%0.6f\n",
            team, rowd[team], N2, t, mflops
        );
    }

    free(rowd);
    free(tGroups);

    double etime = (tend - tstart) / (double)numInvocations;
    double pSize = (double)N1*(double)N2;
    double speed = 5.0 * pSize * log2(pSize) / etime;
    double mflops = speed * 1e-06;

    printf(
        "m=%d, n=%d, time(sec)=%0.6f, Speed(mflops)=%0.6f\n",
        N1, N2, etime, mflops
    );

    exit(EXIT_SUCCESS);
}

/*----------------------------------------------------------------*/
