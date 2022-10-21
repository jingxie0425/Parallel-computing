CXX = charmc
CXXFLAGS = -O2

all: glife2d

glife2d.decl.h : glife2d.ci
	$(CXX) glife2d.ci

glife2d.o : glife2d.C glife2d.decl.h
	$(CXX) -c $(CXXFLAGS) glife2d.C

glife2d : glife2d.o
	$(CXX) -language charm++ -o glife2d glife2d.o

clean:
	rm -f glife2d *.decl.h *.def.h *.o charmrun