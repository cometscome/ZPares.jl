include ../make.inc
include ../make.sufrule

ifdef USE_MPI
ifneq ($(USE_MPI), 0)
DEF_MPI = -DMPI
endif
endif

OBJS = zpares.o zpares_aux.o
OBJ_MUMPS = zpares_mumps.o

ifdef USE_MUMPS
ifneq ($(USE_MUMPS), 0)
ZP_INCLUDE = -I$(MUMPS_DIR)/include/
TARGET_MUMPS = zpares_mumps
endif
endif

all : zpares $(TARGET_MUMPS)

zpares : $(OBJS)
	ar r ../lib/libzpares.a $(OBJS); cp zpares.mod ../include

zpares.o : zpares_aux.o zpares_rci_inc.f90 zpares_rci_sub_inc.f90 zpares_dense_inc.f90

zpares_aux.o : zpares_aux_inc.f90

zpares_mumps : zpares_mumps.o
	ar r ../lib/libzpares_mumps.a $(OBJ_MUMPS); cp zpares_mumps.mod ../include

zpares_mumps.o : zpares_mumps_inc.f90 zpares.o

clean :
	rm -f *.o *.mod ../include/zpares.mod ../include/zpares_mumps.mod ../lib/libzpares.a ../lib/libzpares_mumps.a
