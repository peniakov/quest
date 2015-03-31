include ../../make.inc

FLAGS=$(FC_FLAGS) -I../../SRC

.SUFFIXES : .mod .o .F90

test:   ../../$(DQMCLIB) ggeom.F90
	$(FC) $(FLAGS) -o ggeom ggeom.F90 ../../$(DQMCLIB) $(LIB)

#BetSal:  ../../$(DQMCLIB) ggeom_bethe_salpeter.F90
#	$(FC) $(FLAGS) -o BetSal ggeom_bethe_salpeter.F90 ../../$(DQMCLIB) $(LIB)

clean:
	rm -f ggeom
