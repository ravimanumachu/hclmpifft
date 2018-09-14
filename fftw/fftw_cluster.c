
/*----------------------------------------------------------------*/

#include "hclfftw.h"
#include <omp.h>

/*----------------------------------------------------------------*/

fftw_plan fftw1d_init_plan(
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

    return fftw_plan_many_dft(
               rank, s, howmany,
               X, inembed, istride, idist,
               Y, onembed, ostride, odist,
               sign, FFTW_ESTIMATE);
}



// RICO:

fftw_plan fftw2d_init_plan(
                           const int sign,
                           const int m, const int n,
                           fftw_complex* X, fftw_complex* Y
                           )
{
    /*
     * A 2D row FFTs of size 'm' x 'n'
     */
    int rank = 2;
    int howmany = 1;
    int s[] = {m,n};
    int idist = 0;
    int odist = 0;
    int istride = 1;
    int ostride = 1;
    int *inembed = s;
    int *onembed = s;
    
    return fftw_plan_many_dft(
                              rank, s, howmany,
                              X, inembed, istride, idist,
                              Y, onembed, ostride, odist,
                              sign, FFTW_ESTIMATE);
}


fftw_plan fftw2d_init_plan_basic(
                                 const int sign,
                                 const int m, const int n,
                                 fftw_complex* X, fftw_complex* Y
                                 )
{
    /*
     * A 2D row FFTs of size 'm' x 'n'
     */
    int rank = 2;
    int howmany = 1;
    int s[] = {m,n};
    int idist = 0;
    int odist = 0;
    int istride = 1;
    int ostride = 1;
    int *inembed = s;
    int *onembed = s;
    
    return fftw_plan_dft_2d(
                            m, n, X, Y, FFTW_FORWARD, FFTW_ESTIMATE);
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

    fftw_plan plan = fftw1d_init_plan(
                        sign, *m, n,
                        lMatrix, lMatrix);

    fftw_execute(plan);

    fftw_destroy_plan(plan);

    tGroups[0] += omp_get_wtime() - ts;

    return 0;
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
    fftw_plan plan1 = fftw1d_init_plan(
                        sign, m[0], n,
                        lMatrix, lMatrix);

    fftw_plan plan2 = fftw1d_init_plan(
                        sign, m[1], n,
                        lMatrix, lMatrix);

#pragma omp parallel sections num_threads(nThreadsPerGroup)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan1);
       fftw_destroy_plan(plan1);
       tGroups[0] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan2);
       fftw_destroy_plan(plan2);
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
    fftw_plan plan1 = fftw1d_init_plan(
                        sign, m[0], n,
                        lMatrix, lMatrix);

    fftw_plan plan2 = fftw1d_init_plan(
                        sign, m[1], n,
                        lMatrix, lMatrix);

    fftw_plan plan3 = fftw1d_init_plan(
                        sign, m[2], n,
                        lMatrix, lMatrix);

    fftw_plan plan4 = fftw1d_init_plan(
                        sign, m[3], n,
                        lMatrix, lMatrix);

#pragma omp parallel sections num_threads(nThreadsPerGroup)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan1);
       fftw_destroy_plan(plan1);
       tGroups[0] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan2);
       fftw_destroy_plan(plan2);
       tGroups[1] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan3);
       fftw_destroy_plan(plan3);
       tGroups[2] += omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan4);
       fftw_destroy_plan(plan4);
       tGroups[3] += omp_get_wtime() - ts;
    }

}

    return 0;
}


/*----------------------------------------------------------------*/

int
fftw2dlocal6(
             const int sign,
             const int* m, const int n,
             const int nThreadsPerGroup,
             fftw_complex* lMatrix,
             double* tGroups
             )
{
    fftw_plan plan1 = fftw1d_init_plan(
                                       sign, m[0], n,
                                       lMatrix, lMatrix);
    
    fftw_plan plan2 = fftw1d_init_plan(
                                       sign, m[1], n,
                                       lMatrix, lMatrix);
    
    fftw_plan plan3 = fftw1d_init_plan(
                                       sign, m[2], n,
                                       lMatrix, lMatrix);
    
    fftw_plan plan4 = fftw1d_init_plan(
                                       sign, m[3], n,
                                       lMatrix, lMatrix);

    fftw_plan plan5 = fftw1d_init_plan(
                                       sign, m[4], n,
                                       lMatrix, lMatrix);

    fftw_plan plan6 = /*fftw2d_init_plan_basic*/fftw1d_init_plan(
                                       sign, m[5], n,
                                       lMatrix, lMatrix);
    
#pragma omp parallel sections num_threads(nThreadsPerGroup)
    {
#pragma omp section
        {
            double ts = omp_get_wtime();
            fftw_execute(plan1);
            fftw_destroy_plan(plan1);
            tGroups[0] += omp_get_wtime() - ts;
        }
        
#pragma omp section
        {
            double ts = omp_get_wtime();
            fftw_execute(plan2);
            fftw_destroy_plan(plan2);
            tGroups[1] += omp_get_wtime() - ts;
        }
        
#pragma omp section
        {
            double ts = omp_get_wtime();
            fftw_execute(plan3);
            fftw_destroy_plan(plan3);
            tGroups[2] += omp_get_wtime() - ts;
        }
        
#pragma omp section
        {
            double ts = omp_get_wtime();
            fftw_execute(plan4);
            fftw_destroy_plan(plan4);
            tGroups[3] += omp_get_wtime() - ts;
        }

#pragma omp section
        {
            double ts = omp_get_wtime();
            fftw_execute(plan5);
            fftw_destroy_plan(plan5);
            tGroups[4] += omp_get_wtime() - ts;
        }

#pragma omp section
        {
            double ts = omp_get_wtime();
            fftw_execute(plan6);
            fftw_destroy_plan(plan6);
            tGroups[5] += omp_get_wtime() - ts;
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

//   RICO:

/*----------------------------------------------------------------*/

int fftwlocal_process(
              const int sign,
              const int *m,
              const int n,
              const int nThreadsPerGroup,
              const int nThreadGroups,
              fftw_complex *lMatrix,
              double* tGroups
              )
{

    if (nThreadGroups == 4) {
        
        fftw2dlocal6(sign, m, n, nThreadsPerGroup, lMatrix, tGroups);
        
    } else {
        fprintf(stderr, "ERROR: incorrect number of groups: %d\n", nThreadGroups);
    }
    
    return 0;
}



/*----------------------------------------------------------------*/
