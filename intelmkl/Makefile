
####################################################
#
#    To create mpifftw...
#
####################################################

CC              =       icc
MPICC           =       mpiicc
COPTFLAGS       =       -O3
OMPFLAGS        =       -qopenmp

#---------------------------------------------------#

MKLHOME         =       /opt/intel/mkl
MKLINCLUDES     =       ${MKLHOME}/include
MKLFFTINCLUDES  =       ${MKLHOME}/include/fftw
MKLLIBDIR       =       ${MKLHOME}/lib/intel64
MKLCOMPILERDIR  =       /opt/intel/compilers_and_libraries/linux/lib/intel64
MKLFLAGS        =       ${COPTFLAGS} ${OMPFLAGS} -m64 -Wall -Werror -Wl,--no-as-needed

#---------------------------------------------------#

COMMONINCLUDES	=	../common
EXES		=	mpifftw
SRC		=	mpifftw.c \
			../common/fpm.c \
			../common/mpicomm.c \
			../common/initmatrix.c \
			fftw.c

#---------------------------------------------------#

all: mpifftw

mpifftw:
	${MPICC} ${COPTFLAGS} ${OMPFLAGS} ${MKLFLAGS} \
		-I${COMMONINCLUDES} \
		-I${MKLINCLUDES} -I${MKLFFTINCLUDES} -o $@ ${SRC} \
		-L${MKLLIBDIR} -lmkl_rt \
		-L${MKLCOMPILERDIR} -liomp5 -lpthread -ldl

clean:
	rm -f ${EXES}

#---------------------------------------------------#
