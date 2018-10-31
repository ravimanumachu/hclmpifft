
/*----------------------------------------------------------------*/

#ifndef _HCLFPM_HH
#define _HCLFPM_HH

/*----------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------*/

int determineRowDistributionLB(
    const int n, const int p,
    unsigned int* rowd
);

/*----------------------------------------------------------------*/

int determineRowDistributionLIMB(
    const int n, const int p,
    const unsigned int* npoints,
    const unsigned int* psizes, const double* etimes,
    unsigned int* dOpt
);

/*----------------------------------------------------------------*/

#endif /* _HCLFPM_HH */

/*----------------------------------------------------------------*/
