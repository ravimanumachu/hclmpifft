
#-------------------------------------------------------#

CC		=	gcc
COPTFLAGS       =       -O3 -g -pthread

#---------------------------------------------------#

FFTWHOME        =       ../FFTW_INSTALL
FFTWINCLUDES    =       -I${FFTWHOME}/include \
			-I${WATTSUPHOME}/include -I${GSLHOME}/include
FFTWLIBDIR      =       ${FFTWHOME}/lib
WATTSUPLIBS	=	-L${WATTSUPHOME}/lib -lhclwattsup -L${GSLHOME}/lib -lgsl -lgslcblas

#-------------------------------------------------------#

EXES=fftw

#-------------------------------------------------------#

all: ${EXES}

fftw:
	@echo "--------------------------------------------"
	@echo "FFTW3..."
	@echo "--------------------------------------------"
	@echo
	gcc ${COPTFLAGS} ${FFTWINCLUDES} -c fftw3.c
	g++ ${COPTFLAGS} ${FFTWINCLUDES} \
		-o fftw fftw.cpp fftw3.o \
		-L${FFTWLIBDIR} -lfftw3 -lm \
		${WATTSUPLIBS}

clean:
	rm -f ${EXES} *.o

#-------------------------------------------------------#
