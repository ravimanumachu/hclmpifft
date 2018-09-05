
/*----------------------------------------------------------------*/

#ifndef _HCLFFTCOMM_HH
#define _HCLFFTCOMM_HH

/*----------------------------------------------------------------*/

#include <fftw3.h>

/*----------------------------------------------------------------*/

int mpitranspose(
    const int p, 
    const int *rowd, 
    const int ngroups,
    const int *rowdlocal,
    const int n, 
    fftw_complex* gMatrix
);

/*----------------------------------------------------------------*/

#endif /* _HCLFFTCOMM_HH */

/*----------------------------------------------------------------*/
