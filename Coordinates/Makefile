FCOMP=ifort
INC=
LIB=-L/star/lib

LIBS=-lsla
FLAGS=

all: Coordinates Ephem

Coordinates:

	f2py -c --f90exec=$(FCOMP) $(FLAGS) -m fCoord_tpt fCoord_tpt.f90 $(INC) $(LIB) $(LIBS)

Ephem:

	f2py -c --f90exec=$(FCOMP) $(FLAGS) -m fEphem fEphem.f90 $(INC) $(LIB) $(LIBS)
