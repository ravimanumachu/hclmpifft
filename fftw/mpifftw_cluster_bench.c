
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
// Benchmarking distributed FFTW
// Generate several files with times from each process.

int main(int argc, char **argv) {
    
    int           provided;
    unsigned int  me;
    unsigned int  p;
    unsigned int  inputsIncorrect = 0;
    unsigned int  n;
    unsigned int  ngroups;
    unsigned int  nthreadspergroup;
    unsigned int  verbosity;
    unsigned int  NREPS;
    unsigned int  lb;
    
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int threads_ok = provided >= MPI_THREAD_FUNNELED;
    if (threads_ok) {
        if (verbosity && (me == 0)) {
            fprintf(stdout, "Supporting threads MPI_THREAD_FUNNELED or higher\n"); fflush(stdout);
        }
        threads_ok = fftw_init_threads();
    }
    
    // 1) Read parameters and broadcast them
    if (me == 0) {
        if (argc != 7) {
            fprintf(stderr,
                    "Usage: mpirun -np <p> %s "
                    "              <n> <ngroups> <nthreadspergroup>\n"
                    "              <Load-balanced or Load-imbalanced>\n"
                    "              <verbosity <0|1>\n"
                    "              <num. repetitions (performance) <1>\n"
                    "n is the size of the square signal matrix.\n"
                    "ngroups is the number of groups of threads.\n"
                    "nthreadspergroup is the number of threads per group.\n"
                    "ngroups, nthreadspergroup are for local OpenMP computations\n",
                    argv[0]);
            inputsIncorrect = 1;
        }
    }
    
    int rc = MPI_Bcast(&inputsIncorrect, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "[%d]: Problems broadcasting.\n", me);
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
        NREPS            = atoi(argv[6]);
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
    
    rc = MPI_Bcast(&NREPS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "(%d):Problems broadcasting Num. repetitions\n", me);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    
    
    // 2) Plan with threads
    if (threads_ok) {
        fftw_plan_with_nthreads(nthreadspergroup);
    }
    
    
    // 3) Determine row distribution (homogeneous or using FPMs)
    int* rowd = (int*)malloc(p*sizeof(int));
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
    
    // 4) Determine row -local- distribution (homogeneous or using FPMs)
    // Actually not used in the distributed MPI implementation.
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
        printf("Num groups %d, Threads per group %d, %d: My rowd (%d). Rowd local: ",
               ngroups, nthreadspergroup, me, rowd[me]);
        
        for (int g = 0; g < nthreadspergroup; g++) {
            printf("%d ", rowdlocal[g]);
        }
        printf("\n");
    }
    
    
    // 5) Create local matrix (row[me] x N elements)
    fftw_complex* lMatrix;
    int localnelements = rowd[me] * n;
    lMatrix = (fftw_complex*)fftw_malloc(localnelements * sizeof(fftw_complex));
    if (lMatrix == 0) {
        fprintf(stderr, "Memory allocation failure.\n");
        exit(EXIT_FAILURE);
    }
    
    hclFillSignal2D(rowd[me], n, nthreadspergroup * ngroups, lMatrix);
    
    
    
    // 6) BENCHMARKING
    // 6.a) Create sample vectors
    double *comp_times  = (double *) calloc (NREPS, sizeof(double)); // Computation times
    double *comm_times  = (double *) calloc (NREPS, sizeof(double)); // Communication times
    double *comp_times2 = (double *) calloc (NREPS, sizeof(double)); // Computation times
    double *comm_times2 = (double *) calloc (NREPS, sizeof(double)); // Communication times
    double *total_times = (double *) calloc (NREPS, sizeof(double)); // Total times
    
    // 6.b) Substract the cost of an MPI_Wtime (actually, it does not seem meaningful)
    double tw_start = MPI_Wtime();
    for (int k = 0; k < 1000; k++) {
        MPI_Wtime();
    }
    double wtime_time = (MPI_Wtime() - tw_start) / 1000;
    
    // 6.c) Execution
    for (int k = 0; k < NREPS; k++) {
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        // 8.c.I) FFTW (1)
        double tm_comp_start = MPI_Wtime();
        fftwlocal(FFTW_FORWARD, rowdlocal, n, nthreadspergroup, ngroups, lMatrix, tGroups);
        double tm_comp_end = MPI_Wtime();
        comp_times[k] = tm_comp_end - tm_comp_start;
        
        // 8.c.II) Transpose (1)
        double tm_comm_start = MPI_Wtime();
        mpitranspose(p, rowd, ngroups, nthreadspergroup, rowdlocal, n, lMatrix);
        double tm_comm_end = MPI_Wtime();
        comm_times[k] = tm_comm_end - tm_comm_start;
        
        // 8.c.III) FFTW (2)
        double tm_comp_start2 = MPI_Wtime();
        fftwlocal(FFTW_FORWARD, rowdlocal, n, nthreadspergroup, ngroups, lMatrix, tGroups);
        double tm_comp_end2 = MPI_Wtime();
        comp_times2[k] = tm_comp_end2 - tm_comp_start2;
        
        // 8.c.IV) Transpose (2)
        double tm_comm_start2 = MPI_Wtime();
        mpitranspose(p, rowd, ngroups, nthreadspergroup, rowdlocal, n, lMatrix);
        double tm_comm_end2 = MPI_Wtime();
        comm_times2[k] = tm_comm_end2 - tm_comm_start2;
        
        total_times[k] = (tm_comm_end2 - tm_comp_start) - (8 * wtime_time);
        
        // Tests and progress
        if (verbosity && me == 0) {
            fprintf(stdout, "Rep: %d of %d\n", k+1, NREPS);
            fflush(stdout);
        }
        
    }
    
    
    // 6.d) Save times
    FILE* process_times_file = NULL;
    char name_file[512];
    
    // FORMAT:  process_times_[Rank]_[P]_[M]_[NxN]_[Q]_[thrs].txt
    sprintf(name_file, "process_times_%d_%d_%d_%dx%d_%d_%d.txt", me, p, p/ngroups, n, n, ngroups, nthreadspergroup);
    process_times_file = fopen(name_file, "w");
    
    fprintf(process_times_file, "#NREPS\t %d\n", NREPS);
    fprintf(process_times_file, "#MPI_Wtime\t %0.12f\n", wtime_time);
    
    fprintf(process_times_file, "#NRep\t comp_times(1)\t comm_times(1)\t comp_times(2)\t comm_times(2)\t TOT_times\n");
    
    for (int k = 0; k < NREPS; k++) {
        fprintf(process_times_file, "%d\t %0.12f\t %0.12f\t %0.12f\t %0.12f\t %0.12f\n",
                k,
                comp_times[k],
                comm_times[k],
                comp_times2[k],
                comm_times2[k],
                total_times[k]);
    }
    
    fclose(process_times_file);
    
    // 6.e) Collect data from other processes to show approximate time.
    double t_times[NREPS];
    MPI_Reduce(total_times, &t_times[0], NREPS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) {
        
        int P;
        MPI_Comm_size(MPI_COMM_WORLD, &P);
        for (int i = 0; i < NREPS; i++) {
            t_times[i] = t_times[i] / P;
        }
        
        double time = 0.0;
        for (int i = 0; i < NREPS; i++) {
            time += t_times[i];
        }
        time = time / NREPS;
        fprintf(stdout, "Aprox. Total Time: %0.12f\n", time);
    }
    
    free(comp_times);
    free(comm_times);
    free(comp_times2);
    free(comm_times2);
    free(total_times);

    
    // 7) Cleanup
    if (threads_ok) {
        fftw_cleanup_threads();
    }
    
    free(rowd);
    free(rowdlocal);
    fftw_free(lMatrix);
    free(tGroups);
    
    MPI_Finalize();
    
    exit(EXIT_SUCCESS);
}

