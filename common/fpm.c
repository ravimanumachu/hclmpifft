
/*----------------------------------------------------------------*/

#include "fpm.h"

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
    const size_t* psizes, const double* speeds,
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
