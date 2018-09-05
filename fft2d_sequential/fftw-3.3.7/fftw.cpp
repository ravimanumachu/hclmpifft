
/*-----------------------------------------------------------*/

#include <sys/time.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>
#include <mkl.h>

/*-----------------------------------------------------------*/

#include <wattsup.hpp>
#include <wattsupmacros.hpp>

/*-----------------------------------------------------------*/

#ifdef __cplusplus 
extern "C" {
#endif

int executefftw(const int m, const int n, 
    const unsigned int nt,
    const unsigned int verb,
    double* t);

#ifdef __cplusplus 
}
#endif

/*-----------------------------------------------------------*/

int main(int argc, char** argv)
{
    if (argc != 5)
    {
       std::cerr << "Usage: " << argv[0]
                 << " <M> <N> <nthreads>"
                 << " <verbosity>"
                 << std::endl;
       exit(EXIT_FAILURE);
    }

    MKL_INT M = atoi(argv[1]);
    MKL_INT N = atoi(argv[2]);
    unsigned int nt = atoi(argv[3]);
    unsigned int verbosity = atoi(argv[4]);

    hcl::Precision pIn = {
        100, 0.95, 3600, 0.025
    };
    hcl::Precision pOut;
    hcl::Precision* pOutPtr = &pOut;
    double mean, sd;

    bool precisionAchieved;
    HCL_MEASURE_START(
       pIn, pOutPtr, precisionAchieved);

    double t;
    executefftw(M, N, nt, verbosity, &t);

    HCL_MEASURE_STOP(
       t, 1, pIn, pOutPtr,
       mean, sd);

    double pSize = (double)M*(double)N;
    double speed = 2.5 * pSize * log2(pSize) / mean;
    double mflops = speed * 1e-06;

    printf(
       "m=%d, n=%d, time(sec)=%0.6f, speed(mflops)=%0.6f\n",
       M, N, mean, mflops
    );

    exit(EXIT_SUCCESS);
}

/*-----------------------------------------------------------*/
