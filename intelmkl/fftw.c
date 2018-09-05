
/*----------------------------------------------------------------*/

#include "hclfftw.h"
#include <omp.h>

/*----------------------------------------------------------------*/

int fftw1d(
    const int sign,
    const int m, const int n,
    fftw_complex* X, fftw_complex* Y
)
{
    /*
     * 'm' number of 1D row FFTs of size 'n'
     */
    int rank = 1;
    int howmany = m;
    int s[] = {n};
    int idist = n;
    int odist = n;
    int istride = 1;
    int ostride = 1;
    int *inembed = s;
    int *onembed = s;

    fftw_plan my_plan =  fftw_plan_many_dft(
                           rank, s, howmany,
                           X, inembed, istride, idist,
                           Y, onembed, ostride, odist,
                           sign, FFTW_ESTIMATE);

    fftw_execute(my_plan);

    fftw_destroy_plan(my_plan);

    return 0;
}

/*----------------------------------------------------------------*/

int
fftw2dlocal(
    const int sign,
    const int* m, const int n,
    fftw_complex* lMatrix,
    double* tGroups
)
{
    double ts = omp_get_wtime();
    int rc = fftw1d(
               sign, *m, n,
               lMatrix, lMatrix);
    tGroups[0] += omp_get_wtime() - ts;
    return rc;
}

/*----------------------------------------------------------------*/

int
fftw2dlocal2(
    const int sign,
    const int* m, const int n,
    const int nThreadsPerGroup,
    fftw_complex* lMatrix,
    double* tGroups
)
{
#pragma omp parallel sections num_threads(nThreadsPerGroup)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw1d(
           sign, m[0], n,
           lMatrix, lMatrix);
       tGroups[0] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw1d(
           sign, m[1], n,
           lMatrix, lMatrix);
       tGroups[1] += omp_get_wtime() - ts;
    }
}

    return 0;
}

/*----------------------------------------------------------------*/

int
fftw2dlocal4(
    const int sign,
    const int* m, const int n,
    const int nThreadsPerGroup,
    fftw_complex* lMatrix,
    double* tGroups
)
{
#pragma omp parallel sections num_threads(nThreadsPerGroup)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw1d(
           sign, m[0], n,
           lMatrix, lMatrix);
       tGroups[0] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw1d(
           sign, m[1], n,
           lMatrix, lMatrix);
       tGroups[1] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw1d(
           sign, m[2], n,
           lMatrix, lMatrix);
       tGroups[2] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw1d(
           sign, m[3], n,
           lMatrix, lMatrix);
       tGroups[3] += omp_get_wtime() - ts;
    }
}

    return 0;
}

/*----------------------------------------------------------------*/

int fftwlocal(
    const int sign,
    const int* m,
    const int n,
    const int nThreadsPerGroup,
    const int nThreadGroups,
    fftw_complex *lMatrix,
    double* tGroups
)
{
    if (nThreadGroups == 1)
    {
       fftw2dlocal(
         sign,
         m, n,
         lMatrix, tGroups);
    }

    if (nThreadGroups == 2)
    {
       fftw2dlocal2(
         sign,
         m, n,
         nThreadsPerGroup,
         lMatrix, tGroups);
    }

    if (nThreadGroups == 4)
    {
       fftw2dlocal4(
         sign,
         m, n,
         nThreadsPerGroup,
         lMatrix, tGroups);
    }

    return 0;
}

/*----------------------------------------------------------------*/
