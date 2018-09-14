
/*----------------------------------------------------------------*/

#ifndef _TRANSPOSE_HH
#define _TRANSPOSE_HH

/*----------------------------------------------------------------*/

#include <omp.h>
#include <fftw3.h>

/*----------------------------------------------------------------*/

#define min(a,b) ((a < b) ? a : b)



/*----------------------------------------------------------------*/
/*   Local transpose  (RICO)                                      */
/*----------------------------------------------------------------*/

void hcl_transpose_block_hom (
                              fftw_complex* X,
                              const int start, const int end,
                              const int n,
                              const unsigned int nt,
                              const int block_size,
                              const unsigned int verbosity);

void hcl_transpose_homogeneous (fftw_complex* X,
                                const int start, const int end,
                                const int n,
                                const unsigned int nt,
                                const int block_size,
                                const int *rowd,
                                const unsigned int verbosity);
 
    



/*----------------------------------------------------------------*/

void hcl_local_transpose_scalar_block(
    fftw_complex* X1,
    fftw_complex* X2,
    const int i,
    const int j,
    const int n,
    const int block_size,
    const unsigned int verbosity);

void hcl_local_transpose_block(
    fftw_complex* X,
    const int start, const int end,
    const int n,
    const unsigned int nt,
    const int block_size,
    const unsigned int verbosity);

/*----------------------------------------------------------------*/

#endif /* _TRANSPOSE_HH */

/*----------------------------------------------------------------*/
