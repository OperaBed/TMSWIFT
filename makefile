#CPPLINKER = ${PETSC_DIR}/${PETSC_ARCH}/bin/mpicxx
CPPLINKER = mpicxx
#CPPFLAGS  =  -std=c++0x -lgsl -lgslcblas
CPPFLAGS  = -I${TACC_GSL_INC} -L${TACC_GSL_LIB} -std=c++0x -lgsl -lgslcblas

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

.PHONY: all
all: tmswift analysis plot

tmswift: src/tmswift.o chkopts
	${CPPLINKER} -o tmswift src/tmswift.o  ${CPPFLAGS} ${SLEPC_EPS_LIB}
	${RM} src/tmswift.o

analysis: src/analysis.o chkopts
	${CPPLINKER} -o analysis src/analysis.o  ${CPPFLAGS} ${SLEPC_EPS_LIB}
	${RM} src/analysis.o

#analysis:
#	g++ -o analysis src/analysis.cpp -std=c++0x

plot:
	g++ -o plot_wavefunctions src/plot_wavefunctions.cpp -std=c++0x

.PHONY: clean
clean::
	${RM} tmswift plot_wavefunctions analysis

#to change flags, use CXX_FLAGS in /home/user/petsc-3.6.0/arch-linux2-c-debug/lib/petsc/conf


