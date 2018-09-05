
/*----------------------------------------------------------------*/

#include <fftw3.h>
#include "initmatrix.h"

/*----------------------------------------------------------------*/

int hclFillSignal2D(
    const int m,
    const int n,
    const unsigned int nt,
    fftw_complex* signal
)
{
    int p, q;
#pragma omp parallel for shared(signal) private(p) num_threads(nt)
    for (p = 0; p < m; p++)
    {
        for (q = 0; q < n; q++)
        {
            signal[p*n+q][0] = p*n+q+1;
            signal[p*n+q][1] = p*n+q+1;
        }
    }

    return 0;
}

/*----------------------------------------------------------------*/

int hclPrintSignal2D(
    const int m, const int n,
    const fftw_complex* signal
)
{
    int p, q;
    for (p = 0; p < m; p++)
    {
       for (q = 0; q < n; q++)
       {
           if ((q != 0) && ((q % 4) == 0))
           {
              printf(
                "[%0.2lf + i * %0.2lf]\n",
                signal[p*n+q][0],
                signal[p*n+q][1]
              );
              continue;
           }
           printf(
              "[%0.2lf + i * %0.2lf]\t",
              signal[p*n+q][0],
              signal[p*n+q][1]
           );
       }
       printf("\n");
    }
    printf("\n");

    return 0;
}

/*----------------------------------------------------------------*/
