
/*----------------------------------------------------------------*/

#ifndef _HCLFPM_HH
#define _HCLFPM_HH

/*----------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------*/

int determineRowDistributionLB(
    const int n, const int p,
    int* rowd
);

/*----------------------------------------------------------------*/

int determineRowDistributionLIMB(
    const int n, const int p,
    const size_t* psizes, const double* speeds,
    int* rowd
);

/*----------------------------------------------------------------*/

#endif /* _HCLFPM_HH */

/*----------------------------------------------------------------*/
