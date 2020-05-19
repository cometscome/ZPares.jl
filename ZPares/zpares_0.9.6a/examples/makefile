include ../make.inc
include ../make.sufrule

ifdef USE_MPI
ifneq ($(USE_MPI), 0)
EXEC_MPI = ex2 ex3
DEF_MPI = -DMPI
endif
endif

ifdef USE_MUMPS
ifneq ($(USE_MUMPS), 0)
ifdef USE_MPI
ifneq ($(USE_MPI), 0)
DUMMY_OBJ = blacs_scalapack_dummy.o
DUMMY = $(DUMMY_OBJ)
else
DUMMY = -L$(MUMPS_DIR)/libseq -lmpiseq
endif
else
DUMMY = -L$(MUMPS_DIR)/libseq -lmpiseq
endif
endif
endif

ifdef USE_MUMPS
ifneq ($(USE_MUMPS), 0)
EXEC_MUMPS = ex6
DEF_MUMPS = -DMUMPS
LIB_MUMPS = -lzpares_mumps -L$(MUMPS_DIR)/lib -lcmumps -lzmumps -lmumps_common -lpord $(MUMPS_DEPEND_LIBS)
endif
endif

EXEC = ex1 $(EXEC_MPI) ex4 ex5 $(EXEC_MUMPS) mmzpares

ZP_INCLUDE = -I../include
LINK =  -L../lib $(LIB_MUMPS) -lzpares $(LAPACK) $(BLAS)

all: $(EXEC)

ex1: ex1.o
	$(FC) $(LFFLAG) $^ -o $@ $(LINK)
ex2: ex2.o 
	$(FC) $(LFFLAG) $^ -o $@ $(LINK)
ex3: ex3.o 
	$(FC) $(LFFLAG) $^ -o $@ $(LINK)
ex4: ex4.o 
	$(FC) $(LFFLAG) $^ -o $@ $(LINK)
ex5: ex5.o 
	$(FC) $(LFFLAG) $^ -o $@ $(LINK)
ex6: ex6.o $(DUMMY_OBJ)
	$(FC) $(LFFLAG) ex6.o -o $@ $(LINK) $(DUMMY)

mmzpares: mmzpares.o mmio.o mmzpares_inc.f90 $(DUMMY_OBJ)
	$(FC) $(LFFLAG) mmzpares.o mmio.o -o $@ $(LINK) $(DUMMY)
mmzpares.o:
	$(FC) -c $(FFLAG) $(DEF_MPI) $(DEF_MUMPS) $(ZP_INCLUDE) mmzpares.F90

clean:
	rm -f *.o $(EXEC)
