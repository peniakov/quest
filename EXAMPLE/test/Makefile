include ../../make.inc

.SUFFIXES : .mod .o .F90

test:    ../../$(DQMCLIB) test.F90
	$(FC) $(FLAGS) -I../../SRC -o test  test.F90 ../../$(DQMCLIB) $(LIB) 

clean:
	rm -f test

