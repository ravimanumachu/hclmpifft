
/*----------------------------------------------------------------*/

#ifndef _HCLFFTW_HH
#define _HCLFFTW_HH

/*----------------------------------------------------------------*/

#include <fftw3.h>

/*----------------------------------------------------------------*/

int fftwlocal(
    const int sign,
    const int* m,
    const int n,
    const int nThreadsPerGroup,
    const int nThreadGroups,
    fftw_complex *lMatrix,
    double* tGroups
);

/*----------------------------------------------------------------*/

#endif /* _HCLFFTW_HH */

/*----------------------------------------------------------------*/
