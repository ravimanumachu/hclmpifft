
####################################################
#
#    To create mpifftw...
#
####################################################

MPICC           =       mpicc
COPTFLAGS       =       -O3 -pthread
OMPFLAGS        =       -fopenmp

#---------------------------------------------------#

FFTWINCLUDES    =       ${FFTWHOME}/include
FFTWLIBDIR      =       ${FFTWHOME}/lib

#---------------------------------------------------#

COMMONINCLUDES	=	../common
EXES		=	mpifftw
SRC		=	mpifftw.c \
			../common/initmatrix.c

#---------------------------------------------------#

all: mpifftw

mpifftw:
	${MPICC} ${COPTFLAGS} -I${COMMONINCLUDES} \
		-I${FFTWINCLUDES} -o $@ ${SRC} \
		-L${FFTWLIBDIR} \
		-lfftw3_mpi -lfftw3 -lm

clean:
	rm -f ${EXES}

#---------------------------------------------------#
