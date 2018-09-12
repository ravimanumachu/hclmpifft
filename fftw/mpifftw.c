
/*-----------------------------------------------------------*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

/*-----------------------------------------------------------*/

#include "hclfftw.h"
#include "initmatrix.h"
#include "fpm.h"
#include "mpicomm.h"

/*-----------------------------------------------------------*/

int main(int argc, char **argv)
{
    int provided, me, p;

    MPI_Init_thread(
       &argc, &argv, 
       MPI_THREAD_FUNNELED, 
       &provided);

    int threads_ok = provided >= MPI_THREAD_FUNNELED;
    if (threads_ok)
    {
       threads_ok = fftw_init_threads();
    }

    MPI_Comm_rank(
       MPI_COMM_WORLD, &me);

    MPI_Comm_size(
       MPI_COMM_WORLD, &p);

    unsigned int inputsIncorrect = 0;

    if (me == 0)
    {
       if (argc != 6)
       {
          fprintf(
            stderr,
            "Usage: mpirun -np <p> %s "
            "              <n> <ngroups> <nthreadspergroup>\n"
            "              <Load-balanced or Load-imbalanced>\n"
            "              <verbosity <0|1>\n"
            "n is the size of the square signal matrix.\n"
            "ngroups is the number of groups of threads.\n"
            "nthreadspergroup is the number of threads per group.\n"
            "ngroups, nthreadspergroup are for local OpenMP computations\n",
          argv[0]);
          inputsIncorrect = 1;
       }
    }

    int rc = MPI_Bcast(
               &inputsIncorrect,
            1, MPI_UNSIGNED,
            0, MPI_COMM_WORLD);

    if (rc != MPI_SUCCESS)
    {
       printf("(%d):Problems broadcasting w\n", me);
    }

    if (inputsIncorrect)
    {
       MPI_Finalize();
       exit(EXIT_SUCCESS);
    }

    int n, ngroups, nthreadspergroup, verbosity;
    unsigned int lb;

    if (me == 0)
    {
       n = atoi(argv[1]);
       ngroups = atoi(argv[2]);
       nthreadspergroup = atoi(argv[3]);
       lb = atoi(argv[4]);
       verbosity = atoi(argv[5]);
    }

    rc = MPI_Bcast(
            &n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rc != MPI_SUCCESS)
    {
       fprintf(
          stderr,
          "(%d):Problems broadcasting N\n", me
       );
       MPI_Finalize();
       exit(EXIT_SUCCESS);
    }

    rc = MPI_Bcast(
            &ngroups, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rc != MPI_SUCCESS)
    {
       fprintf(
          stderr,
          "(%d):Problems broadcasting ngroups\n", me
       );
       MPI_Finalize();
       exit(EXIT_SUCCESS);
    }

    rc = MPI_Bcast(
            &nthreadspergroup, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rc != MPI_SUCCESS)
    {
       fprintf(
          stderr,
          "(%d):Problems broadcasting nthreadspergroup\n", me
       );
       MPI_Finalize();
       exit(EXIT_SUCCESS);
    }

    rc = MPI_Bcast(
            &lb, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (rc != MPI_SUCCESS)
    {
       fprintf(
          stderr,
          "(%d):Problems broadcasting lb\n", me
       );
       MPI_Finalize();
       exit(EXIT_SUCCESS);
    }

    rc = MPI_Bcast(
            &verbosity, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rc != MPI_SUCCESS)
    {
       fprintf(
          stderr,
          "(%d):Problems broadcasting verbosity\n", me
       );
       MPI_Finalize();
       exit(EXIT_SUCCESS);
    }

    if (threads_ok) 
    {
       fftw_plan_with_nthreads(nthreadspergroup);
    }

    /*
     * Load the signal matrix on the processor 0
     */
    fftw_complex* gMatrix;
    if (me == 0)
    {
       gMatrix = (fftw_complex*)fftw_malloc(
                      n*n*sizeof(fftw_complex));

       if (gMatrix == NULL)
       {
          fprintf(
            stderr,
            "Memory allocation failure.\n"
          );
          exit(EXIT_FAILURE);
       }

       hclFillSignal2D(
          n, n, 
          nthreadspergroup*ngroups,
          gMatrix);

       if (verbosity)
       {
          printf(
            "Global signal matrix before...\n");
          hclPrintSignal2D(
            n, n, gMatrix);
       }
    }

    unsigned int* rowd = (unsigned int*)malloc(
                         p*sizeof(unsigned int));
    if (me == 0)
    {
       if (lb)
       {
          determineRowDistributionLB(
             n, p,
             rowd);
       }
       else
       {
          /*
           * Determine row distribution from FPMs
           */
          unsigned int* npoints;
          unsigned int* psizes;
          double* etimes;
          determineRowDistributionLIMB(
             n, p,
             npoints, psizes, etimes,
             rowd);
       }
    }

    rc = MPI_Bcast(
            rowd, p, MPI_INT, 0, MPI_COMM_WORLD);

    if (rc != MPI_SUCCESS)
    {
       fprintf(
          stderr,
          "(%d):Problems broadcasting ngroups\n", me
       );
       MPI_Finalize();
       exit(EXIT_SUCCESS);
    }

    unsigned int* rowdlocal = (unsigned int*)malloc(
                              ngroups*sizeof(unsigned int));
    if (lb)
    {
       determineRowDistributionLB(
          rowd[me], ngroups, rowdlocal);
    }
    else
    {
       /*
        * Determine row distribution from FPMs
        */
       unsigned int* npoints;
       unsigned int* psizes;
       double* etimes;
       determineRowDistributionLIMB(
          rowd[me], ngroups,
          npoints, psizes, etimes,
          rowdlocal);
    }

    if (verbosity)
    {
       printf(
          "Num groups %d, Threads per group %d, %d: My rowd %d. Rowd local: ",
          ngroups, nthreadspergroup,
          me, rowd[me]);

       int g;
       for (g = 0; g < ngroups; g++)
       {
           printf("%d ", rowdlocal[g]);
       }
       printf("\n");
    }

    fftw_complex* lMatrix;
    int localnelements = rowd[me]*n;
    lMatrix = (fftw_complex*)fftw_malloc(
                   localnelements*sizeof(fftw_complex));

    if (lMatrix == 0)
    {
       fprintf(
         stderr,
         "Memory allocation failure.\n"
       );
       exit(EXIT_FAILURE);
    }

    /*
     * Scatter the global matrix from processor 0
     * to other processors...
     */
    int* localn;
    int* displs;
    if (me == 0)
    {
       localn = (int*)malloc(p*sizeof(int));
       int proc;
       for (proc = 0; proc < p; proc++)
       {
           localn[proc] = rowd[proc]*n;
       }

       displs = (int*)malloc(p*sizeof(int));
       displs[0] = 0;
       for (proc = 1; proc < p; proc++)
       {
           displs[proc] = displs[proc-1] + localn[proc-1];
       }
    }

    MPI_Scatterv(
       gMatrix, localn, displs,
       MPI_C_DOUBLE_COMPLEX,
       lMatrix, localnelements,
       MPI_C_DOUBLE_COMPLEX,
       0,
       MPI_COMM_WORLD);

    double* tGroups = (double*)calloc(ngroups, sizeof(double));

    /*
     * Compute 'm' number of 1D FFTs
     */
    fftwlocal(
       FFTW_FORWARD, rowdlocal, n,
       nthreadspergroup, ngroups, lMatrix, tGroups);

    /*
     * Use efficient mpi transpose
     */
    mpitranspose(
       p, rowd, ngroups, rowdlocal, n, lMatrix);

    fftwlocal(
       FFTW_FORWARD, rowdlocal, n, 
       nthreadspergroup, ngroups, lMatrix, tGroups);

    mpitranspose(
       p, rowd, ngroups, rowdlocal, n, lMatrix);

    if (threads_ok) 
    { 
       fftw_cleanup_threads();
    }

    MPI_Gatherv(
       lMatrix, localnelements,
       MPI_C_DOUBLE_COMPLEX,
       gMatrix, localn, displs,
       MPI_C_DOUBLE_COMPLEX,
       0,
       MPI_COMM_WORLD);

    if (me == 0)
    {
       if (verbosity)
       {
          printf(
            "Global signal matrix after...\n");
          hclPrintSignal2D(
            n, n, gMatrix);
       }

       free(localn);
       free(displs);
       fftw_free(gMatrix);
    }

    free(rowd);
    free(rowdlocal);
    fftw_free(lMatrix);
    free(tGroups);

    MPI_Finalize();

    exit(EXIT_SUCCESS);
}

/*-----------------------------------------------------------*/
