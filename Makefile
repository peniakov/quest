############################################################################
#
#  Program:         QUEST V1.0
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Modified:        8/20/2008
#
############################################################################

include make.inc

all: libblas liblapack lib example_

libblas:
	(cd BLAS; $(MAKE))

liblapack:
	(cd LAPACK; $(MAKE))

lib:
	(cd SRC; $(MAKE))

example_:
	(cd EXAMPLE; $(MAKE))

clean:
	(cd BLAS; $(MAKE) clean)
	(cd LAPACK; $(MAKE) clean)
	(cd SRC; $(MAKE) clean)
	(cd EXAMPLE; $(MAKE) clean)
	(rm -f $(DQMCLIB))
