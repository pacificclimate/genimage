#CXX=g++-3.4
#CC=gcc-3.4
#CXXFLAGS=-O0 -pg -g -Wall
#CXXFLAGS=-pg -g -Wall
CXXFLAGS=-O3 -funroll-loops 
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg -s
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg
LDFLAGS=-L/usr/local/lib
LDLIBS=-lstdc++ -lgd -lpng -ljpeg -lfreetype -lz
genimage: legend.o genimage.o support.o point.o line.o

%.o: %.cc %.h Makefile

clean: 
	rm *.o genimage
