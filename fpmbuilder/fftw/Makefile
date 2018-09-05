
####################################################
#
#    To create mpifftw...
#
####################################################

CC              =       gcc
COPTFLAGS       =       -O3 -pthread
OMPFLAGS        =       -fopenmp

#---------------------------------------------------#

FFTWINCLUDES    =       ${FFTWHOME}/include
FFTWLIBDIR      =       ${FFTWHOME}/lib

#---------------------------------------------------#

COMMONINCLUDES	=	../../common
EXES		=	constructfpm bestlb
SRC		=	../../common/fpm.c \
			../../common/initmatrix.c \
			../../fftw/fftw.c

#---------------------------------------------------#

all: ${EXES}

${EXES}:
	${CC} ${COPTFLAGS} ${OMPFLAGS} -I${COMMONINCLUDES} \
		-I${FFTWINCLUDES} -o $@ $@.c \
		${SRC} -L${FFTWLIBDIR} \
		-lfftw3_mpi -lfftw3 -lm

clean:
	rm -f ${EXES}

#---------------------------------------------------#
