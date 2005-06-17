CXX=g++
#CXXFLAGS=-O0 -pg -g3 -Wall
#CXXFLAGS=-pg -g -Wall
CXXFLAGS=-march=athlon -O3 -funroll-loops 
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg -s
LDFLAGS=-lstdc++ -lgd -lpng -ljpeg
genimage: legend.o genimage.o support.o point.o line.o

%.o: %.cc %.h Makefile

clean: 
	rm *.o genimage
