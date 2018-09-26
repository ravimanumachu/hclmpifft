
/*----------------------------------------------------------------*/

#include <stdlib.h>
#include <mpi.h>
#include "mpicomm.h"
#include "transpose.h"


/*----------------------------------------------------------------*/



static int HOM_transpose (const int p,
                          const int *rowd,
                          const int ngroups,
                          const int nthreadspergroup,
                          const int *rowdlocal,
                          const int n,
                          fftw_complex* gMatrix);


static int HET_transpose (const int p,
                          const int *rowd,
                          const int ngroups,
                          const int nthreadspergroup,
                          const int *rowdlocal,
                          const int n,
                          fftw_complex* gMatrix);



static int isHomogeneous (const int *rowd, const int P);



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
    
    if (isHomogeneous(rowd, p)) {
        
        return HOM_transpose (p, rowd, ngroups, nthreadspergroup, rowdlocal, N, gMatrix);
        
    } else {
        
        return HET_transpose (p, rowd, ngroups, nthreadspergroup, rowdlocal, N, gMatrix);
    }
}




/*----------------------------------------------------------------*/



/* TBI */
static int isHomogeneous (const int *rowd, const int P) {
    
    int n_rows = rowd[0];
    for (int p = 1; p < P; p++) {
        if (n_rows != rowd[p])
            return 0; // Heterogeneous
    }
    
    return 1;
}






static int HOM_transpose (const int p,
                          const int *rowd,
                          const int ngroups /* NOT USED */,
                          const int nthreadspergroup,
                          const int *rowdlocal /* NOT USED */,
                          const int N,
                          fftw_complex* gMatrix)

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
    int rc = MPI_Alltoall(MPI_IN_PLACE, 1, MPI_TYPE_2D, gMatrix, 1, MPI_TYPE_2D, MPI_COMM_WORLD);
    
    //int counts [4] = {1,1,1,1};
    //int displs [4] = {0,1,2,3};
    // NOT IN PLACE int rc = MPI_Alltoallv(sbuf, counts, displs, MPI_TYPE_2D, rbuf, counts, displs, MPI_TYPE_2D, MPI_COMM_WORLD);
    //int rc = MPI_Alltoallv(MPI_IN_PLACE, counts, displs, MPI_TYPE_2D, gMatrix, counts, displs, MPI_TYPE_2D, MPI_COMM_WORLD);
    
    if (rc != MPI_SUCCESS) {
        printf("[%d] Problems in MPI_Alltoall\n", rank);
    }
    
    MPI_Type_free(&MPI_TYPE_2D);
    
    MPI_Type_free(&MPI_FFTW_COMPLEX);
    
    
    /* Local transpose */
    hcl_transpose_homogeneous (gMatrix, 0, rowd[rank],
                               N, nthreadspergroup, 1, rowd, 0);
    
    return (0);
}






//  RICO:   TBD

static int HET_transpose (const int p,
                          const int *rowd,
                          const int ngroups /* NOT USED */,
                          const int nthreadspergroup,
                          const int *rowdlocal /* NOT USED */,
                          const int N,
                          fftw_complex* gMatrix)

{
    
    int  me;
    int  P;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    
    /*
     if (me == 0) {
     fprintf(stderr, "Error: HETERGENEUOS transpose not fully supported yet (TEMPTATIVE IMPLEMENTATION.");
     }
     */
    
    
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
        
        MPI_Pack((char *)gMatrix + (offset * cplx_ext), rowd[p], MPI_TYPE_COL, buf, tsize, &pos, MPI_COMM_WORLD);
        offset += rowd[p];
    }
    
    
    /* TEMPORAL: imprimir datos *
    sleep(me);
    
    fprintf(stdout, "\n______  %d  _____", me);
    fprintf(stdout, "\nScount:   ");
    for (int p = 0; p < P; p++) {
        fprintf(stdout, "%4d", scounts[p]);
    }
    fprintf(stdout, "\nRcount:   ");
    for (int p = 0; p < P; p++) {
        fprintf(stdout, "%4d", rcounts[p]);
    }
    fprintf(stdout, "\nSdispls:  ");
    for (int p = 0; p < P; p++) {
        fprintf(stdout, "%4d", sdispls[p]);
    }
    fprintf(stdout, "\nRdispls:  ");
    for (int p = 0; p < P; p++) {
        fprintf(stdout, "%4d", rdispls[p]);
    }
    fprintf(stdout, "\nbuf:  ");
    for (int i = 0; i < tsize / cplx_ext; i++) {
        fprintf(stdout, "%2.1f", (double)buf[i*2]);
    }
    fprintf(stdout, "\n");
    
    sleep(P - me);
    */
    
    
    
    MPI_Alltoallw (buf,     scounts, sdispls, stypes,
                   gMatrix, rcounts, rdispls, MPI_TYPE_2D,
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


