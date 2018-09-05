
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

COMMONINCLUDES	=	../../common
EXES		=	constructfpm bestlb
SRC		=	../../common/fpm.c \
			../../common/initmatrix.c \
			../../intelmkl/fftw.c

#---------------------------------------------------#

all: constructfpm bestlb

${EXES}:
	${CC} ${COPTFLAGS} ${OMPFLAGS} ${MKLFLAGS} \
		-I${COMMONINCLUDES} \
		-I${MKLINCLUDES} -I${MKLFFTINCLUDES} -o $@ $@.c ${SRC} \
		-L${MKLLIBDIR} -lmkl_rt \
		-L${MKLCOMPILERDIR} -liomp5 -lpthread -ldl

clean:
	rm -f ${EXES}

#---------------------------------------------------#
