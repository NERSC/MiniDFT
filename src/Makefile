
MKLROOT   = /opt/intel/composer_xe_2015.1.133/mkl

FFTW_INCL = -I/opt/intel/composer_xe_2015.1.133/mkl/include/fftw/
#FFTW_LIBS =  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm ${IPM}
FFTW_LIBS   = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
               $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
               -Wl,--end-group -lpthread -lm ${IPM}

SCALAPACK_LIBS = 

DFLAGS = -D__INTEL

CC = cc 
CFLAGS = -O2

FC = ftn 
FFLAGS = -fpp -O3 -xAVX
#FFLAGS = -fpp -O2 -g

LD = ftn
#LD = hpclink ftn
#LDFLAGS = -g -dynamic

#ifort -openmp defines _OPENMP instead of__OPENMP
DFLAGS    += -D__OPENMP

FFLAGS    += -openmp
LDFLAGS   += -openmp

include Makefile.base
