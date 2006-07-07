CXX=g++-4.0
CC=gcc-4.0
#CXXFLAGS=-O0 -pg -g -Wall
#CXXFLAGS=-pg -g -Wall
CXXFLAGS=-g3 -O2 -march=pentium3 -Wall -Wextra
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg -s
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg
LDFLAGS=-L/usr/local/lib
LDLIBS=-lnetcdf_c++ -lnetcdf -lstdc++ -lgd -lpng -ljpeg -lfreetype -lz
OBJS=legend.o genimage.o support.o point.o line.o ConfigFile.o datamanager.o displayer.o canvas.o range.o
genimage: $(OBJS)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# compile and generate dependency info
%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $*.cc -o $*.o
	$(CXX) -MM $(CXXFLAGS) $*.cc > $*.d

install: genimage
	sudo cp genimage /usr/local/bin/

clean: 
	rm *.o genimage
