CXX=g++-4.0
CC=gcc-4.0
#CXXFLAGS=-O0 -pg -g -Wall
#CXXFLAGS=-pg -g -Wall
CXXFLAGS=-g3 -O2 -march=pentium3 -Wall -Wextra
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg -s
#LDFLAGS=-lstdc++ -lgd -lpng -ljpeg
LDFLAGS=-L/usr/local/lib
LDLIBS=-lnetcdf_c++ -lnetcdf -lstdc++ -lgd -lpng -ljpeg -lfreetype -lz
genimage: legend.o genimage.o support.o point.o line.o ConfigFile.o datamanager.o displayer.o canvas.o range.o

%.o: %.cc %.h Makefile

clean: 
	rm *.o genimage
