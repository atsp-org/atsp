F90 = gfortran
F90FLAGS = -Wall -std=legacy -Wextra -fPIC -Wno-unused -Wno-unused-dummy-argument
F90FLAGS += -g -fcheck=all

all: NONH GENCL

NONH: NONH.o ANG.o COM.o NJGRAF.o RAD6.o
	$(F90) $(F90FLAGS) -o $@ $^

GENCL: GENCL.o
	$(F90) $(F90FLAGS) -o $@ $^

%.o: %.f
	$(F90) $(F90FLAGS) -c $<
