CXX=g++
CC=gcc
#CXXFLAGS=-O0 -g -Wall
CXXFLAGS=-O2 -g -Wall -Wextra -std=c++0x
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg -s
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg
LDFLAGS=-L/usr/local/lib
LDLIBS=-lnetcdf_c++ -lnetcdf -lstdc++ -lgd -lfreetype -lz -lproj -lgdal -lm
OBJS=legend.o genimage.o support.o point.o line.o ConfigFile.o datamanager.o interpdatamanager.o displayer.o canvas.o range.o subset_derived_dm.o derived_vars.o
genimage: $(OBJS)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# compile and generate dependency info
%.o: %.cc Makefile
	$(CXX) -c $(CXXFLAGS) $*.cc -o $*.o
	$(CXX) -MM $(CXXFLAGS) $*.cc > $*.d

install: genimage
	sudo cp genimage /usr/local/bin/

clean: 
	rm -f *.o genimage
