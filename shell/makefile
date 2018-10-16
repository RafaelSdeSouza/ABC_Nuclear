# Define these macros to match with your project's filenames
PROGS = nucleo

FC = gfortran
FORTFLAGS = -O3 -fPIC
#LIBS = -lmkl_rt -liomp5
#	-O3 optimises code
#	-r8 defauls reals to KIND=8 (double precision)
#	-zero zeroes uninitialised variables
#	-CB checks on runtime array bounds
#	"-warn all" prints the warnings

LIBS = -llapack -lblas

all:$(PROGS)
	echo "Done"

nucleo: nucleo.o ma48d.o ddeps.o gearnet.o 
	$(FC) $(FORTFLAGS) $(LIBS) -o $@ $^

%.o: %.f
	$(FC) $(FORTFLAGS) -c $^

nucleo.o: nucleo.f comnuc.inc comgear.inc
	$(FC) $(FORTFLAGS) -c nucleo.f
gearnet.o: gearnet.f comnuc.inc comgear.inc
	$(FC) $(FORTFLAGS) -c gearnet.f

# `make clean' removes extra files that crop up during development
clean:
	-rm *~
	-rm *.o
	-rm $(PROGS)
