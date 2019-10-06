ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX = g++-2.95.3
CINT = $(ROOTSYS)/bin/rootcint
CXXFLAGS = -g $(ROOTCFLAGS)
DEPEND = $(CXX) -MM
DEL = rm -rf

LD = g++-2.95.3
LDFLAGS =  -O4 $(ROOTLIBS)

ROOTHEADERS       := 
SOURCES           := $(wildcard *.cc)
ROOTSRC           := $(ROOTHEADERS:.hh=_root.C)
OBJECTS           := $(SOURCES:.cc=.o) 
ROOTOBJ           := $(ROOTSRC:.C=.o)
EXESRC            := $(wildcard *.cpp)
EXECUTABLES       := $(EXESRC:.cpp=.exe)

all: $(EXECUTABLES)

%.exe : %.cpp $(OBJECTS) $(ROOTOBJ); ${LD} $(CXXFLAGS) $< $(OBJECTS) $(ROOTOBJ) $(LDFLAGS) -o  $@

%.o : %.cc; $(CXX) $(CXXFLAGS) -c $< -o $@

%.o : %.C; $(CXX) $(CXXFLAGS) -c $< -o $@


#%_root.cc : %.hh; $(CINT) -f $@ -c $<

%.d : %.cc; $(DEPEND) $(CXXFLAGS) $< > $@

%.d : %.C; $(DEPEND) $(CXXFLAGS) $< > $@

%.d : %.cpp; $(DEPEND) $(CXXFLAGS) $< > $@; sed s/.o:/.exe:/g $@ > tmp.sed; mv tmp.sed $@ 

clean: ; $(DEL) *.d; $(DEL) *.o; $(DEL) $(ROOTSRC); $(DEL) *.exe; $(DEL) *_root.h; $(DEL) *_root.C; 

include $(SOURCES:.cc=.d)
include $(ROOTSRC:.C=.d)
include $(EXESRC:.cpp=.d)
