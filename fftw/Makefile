
####################################################
#
#    To create mpifftw...
#
####################################################

CC              =       gcc
MPICC           =       mpicxx
COPTFLAGS       =       -O3 -pthread
OMPFLAGS        =       -fopenmp

#---------------------------------------------------#

FFTWINCLUDES    =       ${FFTWHOME}/include
FFTWLIBDIR      =       ${FFTWHOME}/lib

#---------------------------------------------------#

HOPTINCLUDES	=	../hpopt/HPOPT_INSTALL/include
HOPTLIBS    	=	-L../hpopt/HPOPT_INSTALL/lib -lhpopt

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
	${MPICC} ${COPTFLAGS} ${OMPFLAGS} -I${COMMONINCLUDES} \
		-I${HOPTINCLUDES} \
		-I${FFTWINCLUDES} -o $@ \
		${SRC} -L${FFTWLIBDIR} \
		-lfftw3_mpi -lfftw3 \
		${HOPTLIBS} -lm 

clean:
	rm -f ${EXES}

#---------------------------------------------------#
