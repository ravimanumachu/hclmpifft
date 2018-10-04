
/*-----------------------------------------------------------*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

/*-----------------------------------------------------------*/

#include "../common/hclfftw.h"
#include "../common/initmatrix.h"
#include "fpm_cluster.h"
#include "mpicomm.h"

/*-----------------------------------------------------------*/

// A simple program executing a 2D transpose
int main(int argc, char **argv) {
    
    int           provided;
    unsigned int  me;
    unsigned int  p;
    unsigned int  inputsIncorrect = 0;
    unsigned int  n;
    unsigned int  ngroups;
    unsigned int  nthreadspergroup;
    unsigned int  verbosity;
    unsigned int  lb;
    
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    // 1) Evaluate the support for threads in MPI
    int threads_ok = provided >= MPI_THREAD_FUNNELED;
    if (threads_ok) {
        if (verbosity && (me == 0)) {
            fprintf(stdout, "Supporting threads MPI_THREAD_FUNNELED or higher\n"); fflush(stdout);
        }
        threads_ok = fftw_init_threads();
    }
    
    
    // 2) Read parameters and broadcast them
    if (me == 0) {
        if (argc != 6) {
            fprintf(stderr,
                    "Usage: mpirun -np <p> %s "
                    "              <N> <ngroups> <nthreadspergroup>\n"
                    "              <Load-balanced or Load-imbalanced <0|1>\n"
                    "              <verbosity <0|1>\n"
                    "N is the size of the square signal matrix.\n"
                    "ngroups is the number of groups of threads.\n"
                    "nthreadspergroup is the number of threads per group.\n"
                    "ngroups, nthreadspergroup are for local OpenMP computations\n",
                    argv[0]);
            inputsIncorrect = 1;
        }
    }
    
    // 3) Broadcast parameters
    int rc = MPI_Bcast(&inputsIncorrect, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stdout, "[%d]: Problems broadcasting.\n", me);
    }
    if (inputsIncorrect) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if (me == 0) {
        n                = atoi(argv[1]);
        ngroups          = atoi(argv[2]);
        nthreadspergroup = atoi(argv[3]);
        lb               = atoi(argv[4]);
        verbosity        = atoi(argv[5]);
    }
    
    rc = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "(%d):Problems broadcasting N\n", me);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    
    rc = MPI_Bcast(&ngroups, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "(%d):Problems broadcasting ngroups\n", me);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    
    rc = MPI_Bcast(&nthreadspergroup, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "(%d):Problems broadcasting nthreadspergroup\n", me);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    
    rc = MPI_Bcast(&lb, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "(%d):Problems broadcasting lb\n", me);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    
    rc = MPI_Bcast(&verbosity, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "(%d):Problems broadcasting verbosity\n", me);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    
    
    // 4) Plan with threads
    if (threads_ok) {
        fftw_plan_with_nthreads(nthreadspergroup);
    }
    
    
    // 5) Determine row distribution (homogeneous or heterogeneous using FPMs)
    int* rowd = (int *) malloc (p * sizeof(int));
    if (me == 0) {
        if (lb) {
            determineRowDistributionLB(n, p, rowd);
        } else {
            size_t* psizes;
            double* speeds;
            determineRowDistributionLIMB(n, p, psizes, speeds, rowd);
        }
    }
    
    rc = MPI_Bcast(rowd, p, MPI_INT, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "(%d):Problems broadcasting ngroups\n", me);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    
    
    // 6) Determine row -local- distribution (homogeneous or heterogeneous using FPMs)
    //  Actually, this value is not used in MPI implementation (TBD)
    int *rowdlocal = (int *) malloc (nthreadspergroup * sizeof(int));
    if (lb) {
        determineRowDistributionLB(rowd[me], nthreadspergroup, rowdlocal);
    } else {
        size_t* psizes;
        double* speeds;
        determineRowDistributionLIMB(rowd[me], nthreadspergroup, psizes, speeds, rowdlocal);
    }
    
    double* tGroups = (double*)calloc(nthreadspergroup, sizeof(double));

    
    if (verbosity) {
        fprintf(stdout, "Num groups %d, Threads per group %d, %d: My rowd (%d). Rowd local: ",
               ngroups, nthreadspergroup, me, rowd[me]);
        
        for (int g = 0; g < nthreadspergroup; g++) {
            fprintf(stdout, "%d ", rowdlocal[g]);
        }
        fprintf(stdout, "\n");
    }
    
    
    // 7) Create local matrix (row[me] x N elements)
    //    Every process allocates space for its rows.
    fftw_complex* lMatrix;
    int localnelements = rowd[me] * n;
    lMatrix = (fftw_complex*)fftw_malloc(localnelements * sizeof(fftw_complex));
    if (lMatrix == 0) {
        fprintf(stderr, "Memory allocation failure.\n");
        exit(EXIT_FAILURE);
    }
    
    if (verbosity) {
        if (me == 0) fprintf(stdout, "Global signal matrix before...\n");
        hclFillSignal2D(rowd[me], n, nthreadspergroup * ngroups, lMatrix);
        hclPrintSignal2D(rowd[me], n, lMatrix);
    }
        
    
    // 8) Execute distributed FFTW
    double t = MPI_Wtime();
    
    int err;
    
    // 8.a) FFTW (1)
    if (err = fftwlocal(FFTW_FORWARD, rowdlocal, n, nthreadspergroup, ngroups, lMatrix, tGroups)) {
        fprintf(stderr, "ERROR in fftwlocal. \n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }
    
    // 8.b) Transpose (1)
    if (err = mpitranspose(p, rowd, ngroups, nthreadspergroup, rowdlocal, n, lMatrix)) {
        fprintf(stderr, "ERROR in transpose. \n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    // 8.c) FFTW (2)
    if (err = fftwlocal(FFTW_FORWARD, rowdlocal, n, nthreadspergroup, ngroups, lMatrix, tGroups)) {
        fprintf(stderr, "ERROR in fftwlocal. \n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }
    
    // 8.d) Transpose (2)
    if (err = mpitranspose(p, rowd, ngroups, nthreadspergroup, rowdlocal, n, lMatrix)) {
        fprintf(stderr, "ERROR in transpose. \n");
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    t = MPI_Wtime() - t;
    

    if (verbosity) {
        if (me == 0) printf("Global signal matrix after...\n");
        hclPrintSignal2D(rowd[me], n, lMatrix);
    }

    // 9) Cleanup
    if (threads_ok) {
        fftw_cleanup_threads();
    }
    
    free(rowd);
    free(rowdlocal);
    fftw_free(lMatrix);
    free(tGroups);
    
    if (me == 0) {
        fprintf(stdout, "Finishing without errors in %4.6f secs.\n", t);
    }
    
    MPI_Finalize();
    
    exit(EXIT_SUCCESS);
}

