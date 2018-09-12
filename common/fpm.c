
/*----------------------------------------------------------------*/

#include "fpm.h"
#include <hcl_hpopta.h>

/*----------------------------------------------------------------*/

int determineRowDistributionLB(
    const int n, const int p,
    int* rowd
)
{
    /*
     * Simple load-balanced distribution...
     */
    int i, remaining = n%p;
    for (i = 0; i < p; i++)
    {
       rowd[i] = n/p;
    }

    for (i = 0; i < remaining; i++)
    {
       rowd[i] += 1;
    }

    return 0;
}

/*----------------------------------------------------------------*/

int determineRowDistributionLIMB(
    const int n, const int p,
    const unsigned int* npoints,
    const unsigned int* psizes, const double* etimes,
    unsigned int* dOpt
)
{
    double eOpt;
    return hcl_hpopta(
                0, n, p,
                npoints, psizes, etimes,
                dOpt, &eOpt
    );
}

/*----------------------------------------------------------------*/
