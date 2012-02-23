F90 = gfortran
F90FLAGS = -Wall -std=legacy -Wextra -fPIC -Wno-unused -Wno-unused-dummy-argument
F90FLAGS += -g -fcheck=all

PROG = NONH

OBJS = NONH.o ANG.o COM.o NJGRAF.o RAD6.o

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) -o $@ $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.f
	$(F90) $(F90FLAGS) -c $<
