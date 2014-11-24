CXX=g++
CC=gcc
#CXXFLAGS=-O0 -g -Wall
CXXFLAGS=-O2 -g -Wall -Wextra
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg -s
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg
LDFLAGS=-L/usr/local/lib
LDLIBS=-lnetcdf_c++ -lnetcdf -lstdc++
DV_OBJS=common.o compute_derived_vars.o ../genimage/derived_vars.o
CLIM_OBJS=common.o create_climatologies.o
D2M_OBJS=common.o daily2monthly.o

all: compute_derived_vars create_climatologies daily2monthly convert_to_scenarios

compute_derived_vars: $(DV_OBJS)

create_climatologies: $(CLIM_OBJS)

daily2monthly: $(D2M_OBJS)

convert_to_scenarios: common.o

common.o: common.h

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# compile and generate dependency info
%.o: %.cc Makefile
	$(CXX) -c $(CXXFLAGS) $*.cc -o $*.o
	$(CXX) -MM $(CXXFLAGS) $*.cc > $*.d

clean: 
	rm -f $(OBJS) compute_derived_vars

../genimage/derived_vars.o:
	make -C ../genimage/