############################################################################
#
#  Program:         QUEST V1.0
#
#  Module:          make.inc
#
#  Purpose:         System dependent configuration
#
#  Modified date:   10/24/2007
#
#
############################################################################

#
#  Compiler and flags
#

# GNU Fortran compiler
FC = gfortran
FC_FLAGS   = -O3
NOOPT_FLAGS = -O0 

# Intel Fortran compiler
# FC = ifort
# FC_FLAGS    = -O3 -warn -openmp
# NOOPT_FLAGS = -O0 -warn -openmp 

# MAGMA library path
# MAGMAPATH = /home/aetdom/magma_1.1.0
# MAGMALIB = -L$(MAGMAPATH)/lib -lmagma -lmagmablas -lmagma 
# MAGMAINC = -I$(MAGMAPATH)/include

# CUDA installation path
# CUDAPATH = /usr/local/cuda

# CUDA compiler and libraries
# NVCC = $(CUDAPATH)/bin/nvcc
# CU_FLAGS = -O3 -Xptxas -v -m 64 -arch sm_20 $(MAGMAINC)
# CUDALIB = -L$(CUDAPATH)/lib64 -lcublas -lcudart -lcuda
# CUDAINC = -I$(CUDAPATH)/include

# GNU C++ compiler
CXX = gcc
CXX_FLAGS = -Wall -O3 $(CUDAINC) $(MAGMAINC)

# Intel C++ compiler
# CXX = icc
# CXX_FLAGS = -O3 -openmp $(CUDAINC) $(MAGMAINC)

# C++ libraries
CXXLIB = -lstdc++ 

# Included BLAS and LAPACK (not optimized)
LAPACKLIB = ../../liblapack.a ../../libblas.a

# Intel MKL library
# MKLPATH   = /opt/intel/mkl/lib/intel64
# LAPACKLIB = -L$(MKLPATH) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

# BLAS/LAPACK libraries provided by the system, works on most Linux and MacOS X
# LAPACKLIB = -llapack -lblas

DQMCLIB     = libdqmc.a

LIB        = $(CXXLIB) $(LAPACKLIB) $(CUDALIB) $(MAGMALIB) 

#
#  Archiver and its flags
#

ARCH   = ar
ARFLAG = cr
RANLIB = ranlib

#
#  Program dependent flags
#
#    -D_CKB    : use checkboard method for matrix B
#    -D_PREFIX : using prefix matrix multiplication
#    -D_OpenBC : use open boundary condition
#

PRG_FLAGS = -DDQMC_ASQRD # -DDQMC_PROFILE -DDQMC_CUDA -D_CKBP -D_SXX -D_QMC_MPI

FLAGS = $(FC_FLAGS) $(PRG_FLAGS)
