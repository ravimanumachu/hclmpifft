
/*----------------------------------------------------------------*/

#include "hclfftw.h"
#include <sys/time.h>
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

/*----------------------------------------------------------------*/

int
fftw2dlocal(
    const int sign,
    const unsigned int* m, const int n,
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
    const unsigned int* m, const int n,
    const int nThreadsPerGroup,
    fftw_complex* lMatrix,
    double* tGroups
)
{
    double etime, tstart, tend;
    struct timeval start, end;
    gettimeofday(&start, NULL);

    fftw_plan plan1 = fftw1d_init_plan(
                        sign, m[0], n,
                        lMatrix, lMatrix);

    fftw_plan plan2 = fftw1d_init_plan(
                        sign, m[1], n,
                        lMatrix + m[0]*n, 
                        lMatrix + m[0]*n);

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;
    etime = (tend - tstart);

    gettimeofday(&start, NULL);

    printf("g 2, t %d: Plan creation time %f\n", 
          nThreadsPerGroup, etime);

#pragma omp parallel sections num_threads(2)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan1);
       fftw_destroy_plan(plan1);
       tGroups[0] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan2);
       fftw_destroy_plan(plan2);
       tGroups[1] += etime + omp_get_wtime() - ts;
    }
}

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;

    printf("g 2, t %d: Plan execution time %f\n", 
           nThreadsPerGroup, (tend - tstart));

    return 0;
}

/*----------------------------------------------------------------*/

int
fftw2dlocal4(
    const int sign,
    const unsigned int* m, const int n,
    const int nThreadsPerGroup,
    fftw_complex* lMatrix,
    double* tGroups
)
{
    double etime, tstart, tend;
    struct timeval start, end;
    gettimeofday(&start, NULL);

    fftw_plan plan1 = fftw1d_init_plan(
                        sign, m[0], n,
                        lMatrix, lMatrix);

    fftw_plan plan2 = fftw1d_init_plan(
                        sign, m[1], n,
                        lMatrix + m[0]*n, 
                        lMatrix + m[0]*n);

    fftw_plan plan3 = fftw1d_init_plan(
                        sign, m[2], n,
                        lMatrix + (m[0]+m[1])*n, 
                        lMatrix + (m[0]+m[1])*n);

    fftw_plan plan4 = fftw1d_init_plan(
                        sign, m[3], n,
                        lMatrix + (m[0]+m[1]+m[2])*n, 
                        lMatrix + (m[0]+m[1]+m[2])*n);

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;
    etime = (tend - tstart);

    gettimeofday(&start, NULL);

    printf("g 4, t %d: Plan creation time %f\n", 
          nThreadsPerGroup, etime);

#pragma omp parallel sections num_threads(4)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan1);
       fftw_destroy_plan(plan1);
       tGroups[0] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan2);
       fftw_destroy_plan(plan2);
       tGroups[1] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan3);
       fftw_destroy_plan(plan3);
       tGroups[2] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan4);
       fftw_destroy_plan(plan4);
       tGroups[3] += etime + omp_get_wtime() - ts;
    }
}

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;

    printf("g 4, t %d: Plan execution time %f\n", 
           nThreadsPerGroup, (tend - tstart));

    return 0;
}

/*----------------------------------------------------------------*/

int
fftw2dlocal5(
    const int sign,
    const unsigned int* m, const int n,
    const int nThreadsPerGroup,
    fftw_complex* lMatrix,
    double* tGroups
)
{
    double etime, tstart, tend;
    struct timeval start, end;
    gettimeofday(&start, NULL);

    fftw_plan plan1 = fftw1d_init_plan(
                        sign, m[0], n,
                        lMatrix, lMatrix);

    fftw_plan plan2 = fftw1d_init_plan(
                        sign, m[1], n,
                        lMatrix + m[0]*n, 
                        lMatrix + m[0]*n);

    fftw_plan plan3 = fftw1d_init_plan(
                        sign, m[2], n,
                        lMatrix + (m[0]+m[1])*n, 
                        lMatrix + (m[0]+m[1])*n);

    fftw_plan plan4 = fftw1d_init_plan(
                        sign, m[3], n,
                        lMatrix + (m[0]+m[1]+m[2])*n, 
                        lMatrix + (m[0]+m[1]+m[2])*n);

    fftw_plan plan5 = fftw1d_init_plan(
                        sign, m[4], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3])*n, 
                        lMatrix + (m[0]+m[1]+m[2]+m[3])*n);

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;
    etime = (tend - tstart);

    gettimeofday(&start, NULL);

    printf("g 5, t %d: Plan creation time %f\n", 
          nThreadsPerGroup, etime);

#pragma omp parallel sections num_threads(5)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan1);
       fftw_destroy_plan(plan1);
       tGroups[0] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan2);
       fftw_destroy_plan(plan2);
       tGroups[1] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan3);
       fftw_destroy_plan(plan3);
       tGroups[2] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan4);
       fftw_destroy_plan(plan4);
       tGroups[3] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan5);
       fftw_destroy_plan(plan5);
       tGroups[4] += etime + omp_get_wtime() - ts;
    }
}

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;

    printf("g 5, t %d: Plan execution time %f\n", 
           nThreadsPerGroup, (tend - tstart));

    return 0;
}

/*----------------------------------------------------------------*/

int
fftw2dlocal10(
    const int sign,
    const unsigned int* m, const int n,
    const int nThreadsPerGroup,
    fftw_complex* lMatrix,
    double* tGroups
)
{
    double etime, tstart, tend;
    struct timeval start, end;
    gettimeofday(&start, NULL);

    fftw_plan plan1 = fftw1d_init_plan(
                        sign, m[0], n,
                        lMatrix, lMatrix);

    fftw_plan plan2 = fftw1d_init_plan(
                        sign, m[1], n,
                        lMatrix + m[0]*n,
                        lMatrix + m[0]*n);

    fftw_plan plan3 = fftw1d_init_plan(
                        sign, m[2], n,
                        lMatrix + (m[0]+m[1])*n,
                        lMatrix + (m[0]+m[1])*n);

    fftw_plan plan4 = fftw1d_init_plan(
                        sign, m[3], n,
                        lMatrix + (m[0]+m[1]+m[2])*n,
                        lMatrix + (m[0]+m[1]+m[2])*n);

    fftw_plan plan5 = fftw1d_init_plan(
                        sign, m[4], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3])*n);

    fftw_plan plan6 = fftw1d_init_plan(
                        sign, m[5], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4])*n);

    fftw_plan plan7 = fftw1d_init_plan(
                        sign, m[6], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5])*n);

    fftw_plan plan8 = fftw1d_init_plan(
                        sign, m[7], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6])*n);

    fftw_plan plan9 = fftw1d_init_plan(
                        sign, m[8], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7])*n);

    fftw_plan plan10 = fftw1d_init_plan(
                        sign, m[9], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8])*n);

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;
    etime = (tend - tstart);

    gettimeofday(&start, NULL);

    printf("g 10, t %d: Plan creation time %f\n", 
          nThreadsPerGroup, etime);

#pragma omp parallel sections num_threads(10)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan1);
       fftw_destroy_plan(plan1);
       tGroups[0] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan2);
       fftw_destroy_plan(plan2);
       tGroups[1] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan3);
       fftw_destroy_plan(plan3);
       tGroups[2] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan4);
       fftw_destroy_plan(plan4);
       tGroups[3] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan5);
       fftw_destroy_plan(plan5);
       tGroups[4] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan6);
       fftw_destroy_plan(plan6);
       tGroups[5] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan7);
       fftw_destroy_plan(plan7);
       tGroups[6] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan8);
       fftw_destroy_plan(plan8);
       tGroups[7] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan9);
       fftw_destroy_plan(plan9);
       tGroups[8] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan10);
       fftw_destroy_plan(plan10);
       tGroups[9] += etime + omp_get_wtime() - ts;
    }
}

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;

    printf("g 10, t %d: Plan execution time %f\n", 
           nThreadsPerGroup, (tend - tstart));

    return 0;
}

/*----------------------------------------------------------------*/

int
fftw2dlocal20(
    const int sign,
    const unsigned int* m, const int n,
    const int nThreadsPerGroup,
    fftw_complex* lMatrix,
    double* tGroups
)
{
    double etime, tstart, tend;
    struct timeval start, end;
    gettimeofday(&start, NULL);

    fftw_plan plan1 = fftw1d_init_plan(
                        sign, m[0], n,
                        lMatrix, lMatrix);

    fftw_plan plan2 = fftw1d_init_plan(
                        sign, m[1], n,
                        lMatrix + m[0]*n,
                        lMatrix + m[0]*n);

    fftw_plan plan3 = fftw1d_init_plan(
                        sign, m[2], n,
                        lMatrix + (m[0]+m[1])*n,
                        lMatrix + (m[0]+m[1])*n);

    fftw_plan plan4 = fftw1d_init_plan(
                        sign, m[3], n,
                        lMatrix + (m[0]+m[1]+m[2])*n,
                        lMatrix + (m[0]+m[1]+m[2])*n);

    fftw_plan plan5 = fftw1d_init_plan(
                        sign, m[4], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3])*n);

    fftw_plan plan6 = fftw1d_init_plan(
                        sign, m[5], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4])*n);

    fftw_plan plan7 = fftw1d_init_plan(
                        sign, m[6], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5])*n);

    fftw_plan plan8 = fftw1d_init_plan(
                        sign, m[7], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6])*n);

    fftw_plan plan9 = fftw1d_init_plan(
                        sign, m[8], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7])*n);

    fftw_plan plan10 = fftw1d_init_plan(
                        sign, m[9], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8])*n);

    fftw_plan plan11 = fftw1d_init_plan(
                        sign, m[10], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9])*n);

    fftw_plan plan12 = fftw1d_init_plan(
                        sign, m[11], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10])*n);

    fftw_plan plan13 = fftw1d_init_plan(
                        sign, m[12], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11])*n);

    fftw_plan plan14 = fftw1d_init_plan(
                        sign, m[13], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12])*n);

    fftw_plan plan15 = fftw1d_init_plan(
                        sign, m[14], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13])*n);

    fftw_plan plan16 = fftw1d_init_plan(
                        sign, m[15], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14])*n);

    fftw_plan plan17 = fftw1d_init_plan(
                        sign, m[16], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15])*n);

    fftw_plan plan18 = fftw1d_init_plan(
                        sign, m[17], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15]+m[16])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15]+m[16])*n);

    fftw_plan plan19 = fftw1d_init_plan(
                        sign, m[18], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15]+m[16]+m[17])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15]+m[16]+m[17])*n);

    fftw_plan plan20 = fftw1d_init_plan(
                        sign, m[19], n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15]+m[16]+m[17]+m[18])*n,
                        lMatrix + (m[0]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[7]+m[8]+m[9]+m[10]+m[11]+m[12]+m[13]+m[14]+m[15]+m[16]+m[17]+m[18])*n);

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;
    etime = (tend - tstart);

    gettimeofday(&start, NULL);

    printf("g 20, t %d: Plan creation time %f\n", 
          nThreadsPerGroup, etime);

#pragma omp parallel sections num_threads(20)
{
    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan1);
       fftw_destroy_plan(plan1);
       tGroups[0] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan2);
       fftw_destroy_plan(plan2);
       tGroups[1] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan3);
       fftw_destroy_plan(plan3);
       tGroups[2] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan4);
       fftw_destroy_plan(plan4);
       tGroups[3] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan5);
       fftw_destroy_plan(plan5);
       tGroups[4] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan6);
       fftw_destroy_plan(plan6);
       tGroups[5] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan7);
       fftw_destroy_plan(plan7);
       tGroups[6] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan8);
       fftw_destroy_plan(plan8);
       tGroups[7] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan9);
       fftw_destroy_plan(plan9);
       tGroups[8] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan10);
       fftw_destroy_plan(plan10);
       tGroups[9] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan11);
       fftw_destroy_plan(plan11);
       tGroups[10] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan12);
       fftw_destroy_plan(plan12);
       tGroups[11] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan13);
       fftw_destroy_plan(plan13);
       tGroups[12] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan14);
       fftw_destroy_plan(plan14);
       tGroups[13] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan15);
       fftw_destroy_plan(plan15);
       tGroups[14] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan16);
       fftw_destroy_plan(plan16);
       tGroups[15] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan17);
       fftw_destroy_plan(plan17);
       tGroups[16] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan18);
       fftw_destroy_plan(plan18);
       tGroups[17] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan19);
       fftw_destroy_plan(plan19);
       tGroups[18] += etime + omp_get_wtime() - ts;
    }

    #pragma omp section
    {
       double ts = omp_get_wtime();
       fftw_execute(plan20);
       fftw_destroy_plan(plan20);
       tGroups[19] += etime + omp_get_wtime() - ts;
    }
}

    gettimeofday(&end, NULL);
    tstart = start.tv_sec + start.tv_usec/1000000.;
    tend = end.tv_sec + end.tv_usec/1000000.;

    printf("g 20, t %d: Plan execution time %f\n", 
           nThreadsPerGroup, (tend - tstart));

    return 0;
}

/*----------------------------------------------------------------*/

int fftwlocal(
    const int sign,
    const unsigned int* m,
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

    if (nThreadGroups == 5)
    {
       fftw2dlocal5(
         sign,
         m, n,
         nThreadsPerGroup,
         lMatrix, tGroups);
    }

    if (nThreadGroups == 10)
    {
       fftw2dlocal10(
         sign,
         m, n,
         nThreadsPerGroup,
         lMatrix, tGroups);
    }

    if (nThreadGroups == 20)
    {
       fftw2dlocal20(
         sign,
         m, n,
         nThreadsPerGroup,
         lMatrix, tGroups);
    }

    return 0;
}

/*----------------------------------------------------------------*/
