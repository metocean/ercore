FFLAGS=-c -O3 -I/usr/local/intel/include -IPA
LDFLAGS=-O3 -IPA -static -L/usr/local/intel/lib -L/usr/local/lib -lnetcdff -lnetcdf -lhdf5 -lhdf5_hl -lhdf5 -lz
#FFLAGS= -c -g -I./ -I/usr/local/intel/include
#LDFLAGS=-g -static /usr/local/intel/lib/libnetcdf.a

FC=ifort
EXE=ercore

MOD= rnkpar.f90 griddata.f90 shoreline.f90 release.f90 output.f90

.SUFFIXES: .f90 .o

OBJ=    $(MOD:.f90=.o)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: $(EXE)

$(EXE): $(OBJ) ercore.o
	$(FC) -o $@ $(OBJ) $(LDFLAGS) ercore.o

test: $(OBJ) test.o
	$(FC) $(LDFLAGS) -o $@ $(OBJ) test.o

clean:
	rm -f $(OBJ) *.mod $(EXE) ercore.o test.o

pyext:
	f2py --fcompiler=gfortran -c flib.pyf shoreline.f90 plume.f90 slipvel.f90 rnkpar.f90 interp.f90
