
/*----------------------------------------------------------------*/

#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include "mpicomm.h"
#include "transpose.h"

// Testing and debugging transpose
const int test_transpose = 0;
const int MATRIX_PRINT_SIZE = 20;


/*----------------------------------------------------------------*/
// Headers of the transpose algorithms


// decide if a transpose is for a homogeneously distributed matrix or not.
static int isHomogeneous (const int *rowd, const int P);


// Homogeneous transpose.
// Using MPI_Alltoall and threads to transpose blocks locally.
static int HOM_transpose (const int P,
                          const int *rowd,
                          const int ngroups,
                          const int nthreadspergroup,
                          const int *rowdlocal,
                          const int n,
                          fftw_complex* signal);

// Heterogeneous transpose
// Use MPI_Alltoallw, with additional memory to temporary store data
static int HET_transpose (const int P,
                          const int *rowd,
                          const int ngroups,
                          const int nthreadspergroup,
                          const int *rowdlocal,
                          const int n,
                          fftw_complex* signal);


// Heterogeneous transpose
// Use MPI_Send/MPI_Recv and minimizes the additional allocated memory
static int HET_transpose_2 (const int P,
                            const int *rowd,
                            const int ngroups /* NOT USED */,
                            const int nthreadspergroup,
                            const int *rowdlocal /* NOT USED */,
                            const int N,
                            fftw_complex* signal);


/*----------------------------------------------------------------*/

// Auxiliary functions: test and debugging


// Fill matrix "signal" with a simple pattern to being able to test transposition.
int mpi_hclFillSignal2D(const int *m, const int n, const unsigned int nt, fftw_complex* signal) {
    
    int  p, q;
    int  me;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    // Starting at ...
    int  prev = 0;
    for (int i = 0; i < me; i++) {
        prev += m[i];
    }
    prev = prev * n + 1;
    int  my_m = m[me];
    
    // Fill with consecutive numbers
#pragma omp parallel for shared(signal) private(p) num_threads(nt)
    for (p = 0; p < my_m; p++)
    {
        for (q = 0; q < n; q++)
        {
            signal[p*n+q][0] = prev + (p*n+q);
            signal[p*n+q][1] = prev + (p*n+q);
        }
    }
    
    return 0;
}


// Print matrix only if short enough
int mpi_hclPrintSignal2D(const int *rowd, const int n, fftw_complex* signal) {
    
    int  me;
    int  P;
    
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    // Too large to print
    if (n > MATRIX_PRINT_SIZE) return -1;
    
    // As matrices to print are quite short, and this function is intended for testing, we can gather
    //  data in a temporal buffer on rank 0. No problems of memory expected.
    int* localn = NULL;
    int* displs = NULL;
    fftw_complex *tmpbuf = NULL;
    int m = 0;

    if (me == 0) {
        
        for (int p = 0; p < P; p++) {
            m += rowd[p];
        }
        
        int complex_size;
        MPI_Type_size(MPI_C_DOUBLE_COMPLEX, &complex_size);
        
        tmpbuf = (fftw_complex *) malloc (m * n * complex_size);
        if (tmpbuf == NULL) {
            fprintf(stderr, "Memory allocation failure.\n");
            exit(EXIT_FAILURE);
        }
        
        localn = (int *) malloc (P * sizeof(int));
        for (int p = 0; p < P; p++) {
            localn[p] = rowd[p] * n;
        }
        
        displs = (int *) malloc (P * sizeof(int));
        displs[0] = 0;
        for (int p = 1; p < P; p++) {
            displs[p] = displs[p-1] + localn[p-1];
        }
    }
    
    MPI_Gatherv(signal, rowd[me] * n, MPI_C_DOUBLE_COMPLEX,
                tmpbuf, localn, displs, MPI_C_DOUBLE_COMPLEX,
                0, MPI_COMM_WORLD);
    
    if (me == 0) {
        
        for (int p = 0; p < m; p++) {
            for (int q = 0; q < n; q++) {
                fprintf(stdout,
                        "[%4.1lf + %4.1lf]",
                        tmpbuf[p*n+q][0],
                        tmpbuf[p*n+q][1]);
            }
            fprintf(stdout, "\n");
        }
    }
    
    if (me == 0) {
        free(localn);
        free(displs);
        free(tmpbuf);
    }
    
    return 0;
}




// Not deep tests: test if a transpose filled with the "consecutive" pattern has been
//  correctly transposed
int mpi_hclTestSignalTranspose2D (const int *m, const int n, fftw_complex* signal) {
    
    
    int     p, q;
    int     me;
    double  v;
    int     my_error = 0;
    int     P;
    
    
    MPI_Comm_rank (MPI_COMM_WORLD, &me);
    MPI_Comm_size (MPI_COMM_WORLD, &P);
    
    int     errors[P];

    // Previous rows
    int prev = 0;
    for (int i = 0; i < me; i++) {
        prev += m[i];
    }
    prev += 1;
    
    for (p = 0; p < m[me]; p++) {
        for (q = 0; q < n; q++)  {
            v = prev + (q * n + p);
            
            if ((signal[p*n+q][0] != v) || (signal[p*n+q][1] != v)) {
                fprintf(stdout, "[%d]: ERROR in [%d,%d (%d)]:  expected: %4.1f  obtained %4.1f / %4.1f \n", me, p, q, prev, v, signal[p*n+q][0], signal[p*n+q][1]);
                return (my_error = 1);
            }
        }
    }
    
    return (EXIT_SUCCESS);
}




/*----------------------------------------------------------------*/

// Algorithm implementation


/* To be improved (TBI) */
static int isHomogeneous (const int *rowd, const int P) {
    
    int n_rows = rowd[0];
    for (int p = 1; p < P; p++) {
        if (n_rows != rowd[p])
            return 0; // Heterogeneous
    }
    
    return 1;
}




int mpitranspose(const int p,
                 const int *rowd,
                 const int ngroups /* NOT USED */,
                 const int nthreadspergroup,
                 const int *rowdlocal /* NOT USED */,
                 const int N,
                 fftw_complex* signal)
{
    
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    int err = 0;
    
    if (test_transpose) {
        mpi_hclFillSignal2D(rowd, N, nthreadspergroup, signal);
        if (me == 0) fprintf(stdout, "Matrix BEFORE transposing ...\n");
        mpi_hclPrintSignal2D(rowd, N, signal);
    }
    
    
    if (isHomogeneous(rowd, p)) {
        err = HOM_transpose (p, rowd, ngroups, nthreadspergroup, rowdlocal, N, signal);
    } else {
        err = HET_transpose (p, rowd, ngroups, nthreadspergroup, rowdlocal, N, signal);
    }

    
    if (test_transpose) {
        if (me == 0) fprintf(stdout, "Matrix AFTER transposing ...\n");
        mpi_hclPrintSignal2D(rowd, N, signal);
        if (mpi_hclTestSignalTranspose2D(rowd, N, signal)) {
            err = 1;
        }
    }
    
    return err;
    
}




/*----------------------------------------------------------------*/


static int HOM_transpose (const int p,
                          const int *rowd,
                          const int ngroups /* NOT USED */,
                          const int nthreadspergroup,
                          const int *rowdlocal /* NOT USED */,
                          const int N,
                          fftw_complex* signal)
{
    
    MPI_Datatype MPI_C_2D;
    int  rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    /* Complex type */
    int dsize = sizeof(fftw_complex);
    
    MPI_Datatype  MPI_FFTW_COMPLEX;
    MPI_Type_contiguous(dsize, MPI_BYTE, &MPI_FFTW_COMPLEX);
    MPI_Type_commit(&MPI_FFTW_COMPLEX);
    
    
    MPI_Datatype MPI_TYPE_COL, MPI_TYPE_2D;
    
    MPI_Aint cplx_lb, cplx_ext;
    cplx_lb = 0;
    MPI_Type_get_extent(MPI_FFTW_COMPLEX, &cplx_lb, &cplx_ext);
    
    MPI_Type_vector(rowd[rank], rowd[rank], N, MPI_FFTW_COMPLEX, &MPI_TYPE_COL);
    MPI_Type_create_resized(MPI_TYPE_COL, 0, rowd[rank] * cplx_ext, &MPI_TYPE_2D);
    MPI_Type_commit(&MPI_TYPE_2D);
    
    int size;
    MPI_Aint lb, extent;
    lb = 0;
    MPI_Type_size(MPI_TYPE_2D, &size);
    MPI_Type_get_extent(MPI_TYPE_2D, &lb, &extent);
    
    
    // Communication:
    // NOT IN PLACE:  MPI_Alltoall(sbuf, 1, MPI_TYPE_2D, rbuf, 1, MPI_TYPE_2D, MPI_COMM_WORLD);
    int rc = MPI_Alltoall(MPI_IN_PLACE, 1, MPI_TYPE_2D, signal, 1, MPI_TYPE_2D, MPI_COMM_WORLD);
    
    if (rc != MPI_SUCCESS) {
        printf("[%d] Problems in MPI_Alltoall\n", rank);
    }
    
    MPI_Type_free(&MPI_TYPE_2D);
    
    MPI_Type_free(&MPI_FFTW_COMPLEX);
    
    
    /* Local transpose */
    hcl_transpose_homogeneous (signal, 0, rowd[rank],
                               N, nthreadspergroup, 1, rowd, 0);
    
    return (0);
}





static int HET_transpose (const int P,
                          const int *rowd,
                          const int ngroups /* NOT USED */,
                          const int nthreadspergroup,
                          const int *rowdlocal /* NOT USED */,
                          const int N,
                          fftw_complex* signal)

{
    int  me;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    /* Complex type */
    int dsize = sizeof(fftw_complex);
    
    MPI_Datatype  MPI_FFTW_COMPLEX;
    MPI_Type_contiguous(dsize, MPI_BYTE, &MPI_FFTW_COMPLEX);
    MPI_Type_commit(&MPI_FFTW_COMPLEX);
    
    
    /* Column type to pack data */
    MPI_Datatype  MPI_TYPE_COLb, MPI_TYPE_COL;
    
    MPI_Aint cplx_lb, cplx_ext;
    cplx_lb = 0;
    MPI_Type_get_extent(MPI_FFTW_COMPLEX, &cplx_lb, &cplx_ext);
    
    MPI_Type_vector(rowd[me], 1, N, MPI_FFTW_COMPLEX, &MPI_TYPE_COLb);
    MPI_Type_create_resized(MPI_TYPE_COLb, 0, cplx_ext, &MPI_TYPE_COL);
    MPI_Type_commit(&MPI_TYPE_COL);
    
    
    /* 2D types to receive data from each process */
    MPI_Datatype  MPI_TYPE_2Db, MPI_TYPE_2D [P];
    MPI_Datatype  stypes [P];
    int  scounts [P];
    int  rcounts [P];
    int  sdispls [P];
    int  rdispls [P];
    int  tsize = 0;
    int  size;
    
    for (int p = 0; p < P; p++) {
        
        /* Types to receive */
        MPI_Type_vector(rowd[me], rowd[p], N, MPI_FFTW_COMPLEX, &MPI_TYPE_2Db);
        MPI_Type_create_resized(MPI_TYPE_2Db, 0, cplx_ext, &MPI_TYPE_2D[p]);
        MPI_Type_commit(&MPI_TYPE_2D[p]);
        
        /* Types to send */
        stypes[p] = MPI_PACKED;
        
        /* Size of the contiguous buffer to send */
        MPI_Pack_size(rowd[p], MPI_TYPE_COL, MPI_COMM_WORLD, &size);
        tsize += size;
        rcounts[p] = 1;
        scounts[p] = size;
        
    }
    
    char *buf = NULL;
    
    if ((buf = (char *) malloc (tsize)) == NULL) {
        fprintf(stderr, "ERROR: allocating memory for temporary buffers\n");
        MPI_Abort(MPI_COMM_WORLD, (-1));
    }
    
    int  pos = 0;
    int  offset = 0;
    
    for (int p = 0; p < P; p++) {
        sdispls[p] = pos;
        rdispls[p] = (p==0) ? 0 : (rowd[p-1] * cplx_ext + rdispls[p-1]);
        
        MPI_Pack((char *)signal + (offset * cplx_ext), rowd[p], MPI_TYPE_COL, buf, tsize, &pos, MPI_COMM_WORLD);
        offset += rowd[p];
    }
    
    
    MPI_Alltoallw (buf,     scounts, sdispls, stypes,
                   signal, rcounts, rdispls, MPI_TYPE_2D,
                   MPI_COMM_WORLD);
    
    free(buf);
    
    MPI_Type_free(&MPI_FFTW_COMPLEX);
    MPI_Type_free(&MPI_TYPE_COLb);
    MPI_Type_free(&MPI_TYPE_COL);
    
    for (int p = 0; p < P; p++) {
        MPI_Type_free(&MPI_TYPE_2D[p]);
    }
    
    return (EXIT_SUCCESS);
}




static int HET_transpose_2 (const int P,
                            const int *rowd,
                            const int ngroups /* NOT USED */,
                            const int nthreadspergroup,
                            const int *rowdlocal /* NOT USED */,
                            const int N,
                            fftw_complex* signal)

{
    int  me;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    /* Complex type */
    int dsize = sizeof(fftw_complex);
    
    MPI_Datatype  MPI_FFTW_COMPLEX;
    MPI_Type_contiguous(dsize, MPI_BYTE, &MPI_FFTW_COMPLEX);
    MPI_Type_commit(&MPI_FFTW_COMPLEX);
    
    
    /* Column type to pack data */
    MPI_Datatype  MPI_TYPE_COLb, MPI_TYPE_COL;
    
    MPI_Aint cplx_lb, cplx_ext;
    cplx_lb = 0;
    MPI_Type_get_extent(MPI_FFTW_COMPLEX, &cplx_lb, &cplx_ext);
    
    MPI_Type_vector(rowd[me], 1, N, MPI_FFTW_COMPLEX, &MPI_TYPE_COLb);
    MPI_Type_create_resized(MPI_TYPE_COLb, 0, cplx_ext, &MPI_TYPE_COL);
    MPI_Type_commit(&MPI_TYPE_COL);
    
    
    /* 2D types to receive data from each process */
    int  scounts [P];
    int  max_size = 0;
    int  size;
    
    for (int p = 0; p < P; p++) {
        
        /* Size of the contiguous buffer to send */
        MPI_Pack_size(rowd[p], MPI_TYPE_COL, MPI_COMM_WORLD, &size);
        max_size = (size > max_size) ? size : max_size;
        scounts[p] = size;
    }
    
    char *buf = NULL;
    
    if ((buf = (char *) calloc (max_size, 1)) == NULL) {
        fprintf(stderr, "ERROR: allocating memory for temporary buffers\n");
        MPI_Abort(MPI_COMM_WORLD, (-1));
    }
    
    int  pos = 0;
    int  offset = 0;
    
    int  MY_ALLTOALL_TAG = 100;
    
    
    for (int p = 0; p < P; p++) {
        
        /* Type to receive */
        MPI_Datatype  MPI_TYPE_2Db, MPI_TYPE_2D;
        
        MPI_Type_vector(rowd[me], rowd[p], N, MPI_FFTW_COMPLEX, &MPI_TYPE_2Db);
        MPI_Type_create_resized(MPI_TYPE_2Db, 0, cplx_ext, &MPI_TYPE_2D);
        MPI_Type_commit(&MPI_TYPE_2D);
        
        pos = 0;
        MPI_Pack((char *)signal + (offset * cplx_ext), rowd[p], MPI_TYPE_COL, buf, max_size, &pos, MPI_COMM_WORLD);
        // TBI: Se pueden crear otro tipo 2D por columnas para para PACK, aunque no creo que mejore.
        //MPI_Pack((char *)signal + (offset * cplx_ext), 1, MPI_TYPE_2D_COL, buf, max_size, &pos, MPI_COMM_WORLD);
        
        if (p == me) {
            
            pos = 0;
            MPI_Unpack(buf, max_size, &pos, (char *)signal + (offset * cplx_ext), 1, MPI_TYPE_2D, MPI_COMM_WORLD);
            
        } else {
            
            MPI_Request req;
            
            MPI_Irecv((char *)signal + (offset * cplx_ext), 1, MPI_TYPE_2D, p, MY_ALLTOALL_TAG, MPI_COMM_WORLD, &req);
            MPI_Send(buf, scounts[p], MPI_PACKED, p, MY_ALLTOALL_TAG, MPI_COMM_WORLD);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
            
        }
        
        offset += rowd[p];
        
        MPI_Type_free(&MPI_TYPE_2Db);
        MPI_Type_free(&MPI_TYPE_2D);
    }
    
    free(buf);
    
    MPI_Type_free(&MPI_FFTW_COMPLEX);
    MPI_Type_free(&MPI_TYPE_COLb);
    MPI_Type_free(&MPI_TYPE_COL);
    
    return (EXIT_SUCCESS);
}


