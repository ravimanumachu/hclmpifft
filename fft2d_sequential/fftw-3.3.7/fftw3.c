
/*-----------------------------------------------------------*/

#include <sys/time.h>
#include "fftw3.h"

/*----------------------------------------------------------------*/

int fillSignal2D(
    const int m,
    const int n,
    fftw_complex* signal
)
{
    int p;
    for (p = 0; p < (m*n); p++)
    {
        signal[p][0] = p+1;
        signal[p][1] = p+1;
    }
 
    return 0;
}

/*-----------------------------------------------------------*/

static void print(fftw_complex *x, int N1, int N2)
{
    int n1, n2, index = 0;

    for (n1 = 0; n1 < N1; n1++)
    {
        for (n2 = 0; n2 < N2; n2++)
        {
            printf("(%lf, %lf) ", x[index][0], x[index][1]);
            index++;
        }
        printf("\n");
    }
}

/*-----------------------------------------------------------*/

int executefftw(const int m, const int n,
    const unsigned int nt,
    const unsigned int verbosity,
    double* t)
{
    fftw_plan forward_plan = 0, backward_plan = 0;

    fftw_complex *x = 0;

    fftw_init_threads();
    fftw_plan_with_nthreads(nt);

    /* Sizes of 2D transform */
    int N[2] = {m, n};
    int H[2] = {-1, -2};

    x  = fftw_malloc(sizeof(fftw_complex)*N[0]*N[1]);

    fillSignal2D(m, n, x);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    forward_plan = fftw_plan_dft(2, N, x, x, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(forward_plan);
    if (verbosity)
    {
       print(x, m, n);
    }

    fftw_destroy_plan(forward_plan);

    backward_plan = fftw_plan_dft(2, N, x, x, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(backward_plan);

    fftw_destroy_plan(backward_plan);

    gettimeofday(&end, NULL);

    double tstart = start.tv_sec + start.tv_usec/1000000.;
    double tend = end.tv_sec + end.tv_usec/1000000.;
    *t = (tend - tstart);

    fftw_free(x);

    fftw_cleanup_threads();

    return 0;
}

/*-----------------------------------------------------------*/
