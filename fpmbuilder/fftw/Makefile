
####################################################
#
#    To create mpifftw...
#
####################################################

CC              =       g++
COPTFLAGS       =       -O3 -pthread
OMPFLAGS        =       -fopenmp

#---------------------------------------------------#

FFTWHOME        =       /lusitania_homes/UniversidadDeExtremadura/juanantonio.rico/hclmpifft-master/FFTW_INSTALL
FFTWINCLUDES    =       ${FFTWHOME}/include
FFTWLIBDIR      =       ${FFTWHOME}/lib

HPOPTAHOME	=	/lusitania_homes/UniversidadDeExtremadura/juanantonio.rico/hclmpifft-master/hpopt/HPOPT_INSTALL
HPOPTAINCLUDES  =	${HPOPTAHOME}/include
HPOPTALIBDIR    =	${HPOPTAHOME}/lib

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
		-I${HPOPTAINCLUDES} \
		-I${FFTWINCLUDES} -o $@ $@.c \
		${SRC} -L${HPOPTALIBDIR} \
		-lhpopt -L${FFTWLIBDIR} -lfftw3 -lm

clean:
	rm -f ${EXES}

#---------------------------------------------------#
