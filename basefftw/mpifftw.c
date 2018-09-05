
/*-----------------------------------------------------------*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

/*-----------------------------------------------------------*/

#include <fftw3-mpi.h>
#include "initmatrix.h"

/*-----------------------------------------------------------*/

int main(int argc, char **argv)
{
    int me, p;

    MPI_Init(&argc, &argv);

    fftw_mpi_init();

    MPI_Comm_rank(
       MPI_COMM_WORLD, &me);

    MPI_Comm_size(
       MPI_COMM_WORLD, &p);

    unsigned int inputsIncorrect = 0;

    if (me == 0)
    {
       if (argc != 3)
       {
          fprintf(
            stderr,
            "Usage: mpirun -np <p> %s <n> <verbosity (0|1)\n"
            "n is the size of the square signal matrix.\n",
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

    int n, verbosity;

    if (me == 0)
    {
       n = atoi(argv[1]);
       verbosity = atoi(argv[2]);
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
          1, gMatrix);

       if (verbosity)
       {
          printf(
            "Global signal matrix before...\n");
          hclPrintSignal2D(
            n, n, gMatrix);
       }
    }

    ptrdiff_t alloc_local, localn, localstart;

    alloc_local = fftw_mpi_local_size_2d(
                      n, n, MPI_COMM_WORLD,
                      &localn, &localstart);

    if (verbosity)
    {
       printf(
         "Me %d: My local rows %zu, alloc local %zu.\n",
         me, localn, alloc_local);
    }

    fftw_complex *lMatrix = fftw_alloc_complex(alloc_local);

    fftw_plan plan = fftw_mpi_plan_dft_2d(
                         n, n, lMatrix, lMatrix, 
                         MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);

    /*
     * Broadcast rows from process 0 to other processes...
     */
    int* rowd;
    if (me == 0)
    {
       rowd = (int*)malloc(p*sizeof(int));
    }

    int localnelements = localn*n;
    MPI_Gather(
       &localnelements, 1,
       MPI_INT,
       rowd, 1,
       MPI_INT, 0,
       MPI_COMM_WORLD);

    int* displs;
    if (me == 0)
    {
       if (verbosity)
       {
          printf(
            "Me %d: Row distribution: ",
            me);
          int proc;
          for (proc = 0; proc < p; proc++)
          {
              printf("%d ", rowd[proc]);
          }
          printf("\n");
       }

       displs = (int*)malloc(p*sizeof(int));
       int proc;
       displs[0] = 0;
       for (proc = 1; proc < p; proc++)
       {
           displs[proc] = displs[proc-1] + rowd[proc-1];
       }
    }

    MPI_Scatterv(
       gMatrix, rowd, displs,
       MPI_C_DOUBLE_COMPLEX,
       lMatrix, localnelements,
       MPI_C_DOUBLE_COMPLEX,
       0,
       MPI_COMM_WORLD);

    if (verbosity)
    {
       printf(
         "Me %d: Signal local matrix before...\n", me);
       hclPrintSignal2D(
         localn, n, gMatrix);
    }

    fftw_execute(plan);

    fftw_destroy_plan(plan);

    MPI_Gatherv(
       lMatrix, localnelements,
       MPI_C_DOUBLE_COMPLEX,
       gMatrix, rowd, displs,
       MPI_C_DOUBLE_COMPLEX,
       0,
       MPI_COMM_WORLD);

    fftw_free(lMatrix);

    if (me == 0)
    {
       if (1)
       {
          printf(
            "Global signal matrix after...\n");
          hclPrintSignal2D(
            n, n, gMatrix);
       }

       free(rowd);
       free(displs);
    }

    if (me == 0)
    {
       fftw_free(gMatrix);
    }

    MPI_Finalize();

    exit(EXIT_SUCCESS);
}

/*-----------------------------------------------------------*/
