FCOMP=ifort
INC=
LIB=-L/star/lib

LIBS=-lsla
FLAGS=

all: Coordinates

Coordinates:

	f2py -c --f90exec=$(FCOMP) $(FLAGS) -m fCoord_tpt fCoord_tpt.f90 $(INC) $(LIB) $(LIBS)

