OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o particleutil.o rotation.o browndyn.o
FC = gfortran
FLAGS = -O3
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -fcheck=all
LIBS = -lblas -llapack
OSGLIBS =  -lblas -llapack -lgfortran -lm
OSGFLAGS = -static 

# executable name
EXE = ../spheretouchBD.exe
TESTEXE = ../testrun.exe
OSGEXE = ../spheretouchBD_osg.exe

.SUFFIXES: .o .f90

# general suffix rules
%.o: %.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: tools/%.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

osg: $(OBJS) main.o
	$(FC) $(FLAGS) $(OSGFLAGS) $(OBJS) main.o -o $(OSGEXE) $(OSGLIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

static: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) -static $(OBJS) main.o -o $(EXE) $(LIBS)

clean:
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o 
particleutil.o: keys.o genutil.o rotation.o
rotation.o: genutil.o
browndyn.o: particleutil.o mt19937.o rotation.o keys.o
