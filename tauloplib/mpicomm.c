
/*----------------------------------------------------------------*/

#include <stdlib.h>
#include <mpi.h>
#include "mpicomm.h"
#include "transpose.h"


/*----------------------------------------------------------------*/



static int HOM_transpose (const int P,
                          const int *rowd,
                          const int ngroups,
                          const int nthreadspergroup,
                          const int *rowdlocal,
                          const int n,
                          fftw_complex* gMatrix);


static int HET_transpose (const int P,
                          const int *rowd,
                          const int ngroups,
                          const int nthreadspergroup,
                          const int *rowdlocal,
                          const int n,
                          fftw_complex* gMatrix);



static int isHomogeneous (const int *rowd, const int P);




static int HET_transpose_2 (const int P,
                            const int *rowd,
                            const int ngroups /* NOT USED */,
                            const int nthreadspergroup,
                            const int *rowdlocal /* NOT USED */,
                            const int N,
                            fftw_complex* gMatrix);

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





static int HET_transpose (const int P,
                            const int *rowd,
                            const int ngroups /* NOT USED */,
                            const int nthreadspergroup,
                            const int *rowdlocal /* NOT USED */,
                            const int N,
                            fftw_complex* gMatrix)

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



static int HET_transpose_2 (const int P,
                          const int *rowd,
                          const int ngroups /* NOT USED */,
                          const int nthreadspergroup,
                          const int *rowdlocal /* NOT USED */,
                          const int N,
                          fftw_complex* gMatrix)

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
        MPI_Pack((char *)gMatrix + (offset * cplx_ext), rowd[p], MPI_TYPE_COL, buf, max_size, &pos, MPI_COMM_WORLD);
        // TBI: Se pueden crear otro tipo 2D por columnas para para PACK, aunque no creo que mejore.
        //MPI_Pack((char *)gMatrix + (offset * cplx_ext), 1, MPI_TYPE_2D_COL, buf, max_size, &pos, MPI_COMM_WORLD);
        
        if (p == me) {
            
            pos = 0;
            MPI_Unpack(buf, max_size, &pos, (char *)gMatrix + (offset * cplx_ext), 1, MPI_TYPE_2D, MPI_COMM_WORLD);
            
        } else {
            
            MPI_Request req;
            
            MPI_Irecv((char *)gMatrix + (offset * cplx_ext), 1, MPI_TYPE_2D, p, MY_ALLTOALL_TAG, MPI_COMM_WORLD, &req);
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


