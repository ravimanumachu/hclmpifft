
/*----------------------------------------------------------------*/

#include <mpi.h>
#include "mpicomm.h"
#include "transpose.h"

/*----------------------------------------------------------------*/

int mpitranspose(
    const int p,
    const int *rowd,
    const int ngroups /* NOT USED */,
    const int nthreadspergroup,
    const int *rowdlocal /* NOT USED */,
    const int N,
    fftw_complex* gMatrix
)
{
    MPI_Datatype MPI_C_2D;
    int  rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //fprintf(stdout, "[%d] Algorithm 3\n", rank); fflush(stdout);
    
    
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
    //MPI_Alltoall(sbuf, 1, MPI_TYPE_2D, rbuf, 1, MPI_TYPE_2D, MPI_COMM_WORLD);
    MPI_Alltoall(MPI_IN_PLACE, 1, MPI_TYPE_2D, gMatrix, 1, MPI_TYPE_2D, MPI_COMM_WORLD);
    //int counts [4] = {1,1,1,1};
    //int displs [4] = {0,1,2,3};
    // ALSO: MPI_Alltoallv(sbuf, counts, displs, MPI_TYPE_2D, rbuf, counts, displs, MPI_TYPE_2D, MPI_COMM_WORLD);
    
    MPI_Type_free(&MPI_TYPE_2D);
    
    MPI_Type_free(&MPI_FFTW_COMPLEX);
    
    
    /* Local transpose */
    hcl_transpose_homogeneous (gMatrix, 0, rowd[rank],
                               N, nthreadspergroup, 1, rowd, 0);
    
    return (0);
}

/*----------------------------------------------------------------*/



//  RICO:   TBD

int hettranspose(
                 const int p,
                 const int *rowd,
                 const int ngroups /* NOT USED */,
                 const int nthreadspergroup,
                 const int *rowdlocal /* NOT USED */,
                 const int N,
                 fftw_complex* gMatrix
                 )
{
    MPI_Datatype MPI_C_2D;
    int  rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //fprintf(stdout, "[%d] Algorithm 3\n", rank); fflush(stdout);
    
    
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
    //MPI_Alltoall(sbuf, 1, MPI_TYPE_2D, rbuf, 1, MPI_TYPE_2D, MPI_COMM_WORLD);
    MPI_Alltoall(MPI_IN_PLACE, 1, MPI_TYPE_2D, gMatrix, 1, MPI_TYPE_2D, MPI_COMM_WORLD);
    //int counts [4] = {1,1,1,1};
    //int displs [4] = {0,1,2,3};
    // ALSO: MPI_Alltoallv(sbuf, counts, displs, MPI_TYPE_2D, rbuf, counts, displs, MPI_TYPE_2D, MPI_COMM_WORLD);
    
    MPI_Type_free(&MPI_TYPE_2D);
    
    MPI_Type_free(&MPI_FFTW_COMPLEX);
    
    
    /* Local transpose */
    hcl_transpose_homogeneous (gMatrix, 0, rowd[rank],
                               N, nthreadspergroup, 1, rowd, 0);
    
    return (0);
}


