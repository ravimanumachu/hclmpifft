
/*----------------------------------------------------------------*/

#ifndef _INITMATRIX_HH
#define _INITMATRIX_HH

/*----------------------------------------------------------------*/

#include <fftw3.h>

/*----------------------------------------------------------------*/

int hclFillSignal2D(
    const int m,
    const int n,
    const unsigned int nt,
    fftw_complex* signal
);

/*----------------------------------------------------------------*/

int hclPrintSignal2D(
    const int m, const int n,
    const fftw_complex* signal
);

/*----------------------------------------------------------------*/

#endif /* _INITMATRIX_HH */

/*----------------------------------------------------------------*/
