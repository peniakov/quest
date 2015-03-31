include ../../make.inc

FLAGS=$(FC_FLAGS) $(PRG_FLAGS) -I../../SRC

.SUFFIXES : .mod .o .F90

mpi:    ../../$(DQMCLIB) verify.F90
	$(FC) $(FLAGS) -o verify verify.F90 ../../$(DQMCLIB) $(LIB) 

clean:
	rm -f verify

