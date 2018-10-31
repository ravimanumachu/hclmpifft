
/*----------------------------------------------------------------*/

#ifndef _HCLFFTCOMM_HH
#define _HCLFFTCOMM_HH

/*----------------------------------------------------------------*/

#include <fftw3.h>

/*----------------------------------------------------------------*/

int mpitranspose(
    const int p, 
    const unsigned int *rowd, 
    const int ngroups,
    const unsigned int *rowdlocal,
    const int n, 
    fftw_complex* gMatrix
);

/*----------------------------------------------------------------*/

#endif /* _HCLFFTCOMM_HH */

/*----------------------------------------------------------------*/
