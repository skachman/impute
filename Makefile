#############################################################################
# Makefile for building: GenSel
# on Darwin OS X
#############################################################################

####### Compiler, tools and options
#COMPILER=INTEL
COMPILER=GCC
#COMPILER=CLANG

DEFINES       = 

MATVECLOC     = matvec

ifeq ($(COMPILER),INTEL)
CXXFLAGS      = -pipe -O3  -g  -arch x86_64 $(DEFINES) #
CC            = icc -fopnemp
CXX           = icpc -fopenmp -std=c++11
LINK          = icpc   -fopenmp 
endif
ifeq ($(COMPILER),GCC)
CC            = gcc -fopenmp
CXX           = g++ -fopenmp -std=c++11
CXXFLAGS      = -pipe -O3 -funsafe-math-optimizations -g  -arch x86_64 $(DEFINES)
LINK          = g++   -fopenmp
endif
ifeq ($(COMPILER),CLANG)
CC   = clang 
CXX  = clang++
LINK          = clang++
CXXFLAGS = -pipe -O3  -g  -arch x86_64 $(DEFINES) 
endif


CFLAGS        = -pipe -g -Wall -W $(DEFINES)

INCPATH       =  -I/opt/local/include/eigen3/ -I/opt/local/include -I. -Iinclude

LFLAGS        = -prebind
LIBS          = $(SUBLIBS)  #/opt/local/lib/libmatvec.a 
AR            = ar cq
RANLIB        = ranlib -s
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = cp -f
COPY_DIR      = cp -f -r
INSTALL_FILE  = $(COPY_FILE)
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = $(COPY_FILE)
SRC_DIR       = src/
HEADER_DIR    = include/
DEL_FILE      = rm -f
SYMLINK       = ln -sf
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   =  obj/

BINARY_DIR    = ~/bin/

####### Files

SOURCES       =  $(SRC_DIR)BayesIM.cpp $(SRC_DIR)BayesIMmt.cpp  $(SRC_DIR)PullRegions.cpp $(SRC_DIR)Configuration.cpp $(SRC_DIR)utility.cpp 


OBJECTS       = $(OBJECTS_DIR)BayesIM.o $(OBJECTS_DIR)BayesIMmt.o $(OBJECTS_DIR)PullRegions.o $(OBJECTS_DIR)Configuration.o $(OBJECTS_DIR)utility.o	

DESTDIR       = 

TARGET        = #bin/impute 

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET) bin/BayesIM bin/BayesIMmt bin/PullRegions # bin/imputeMCMC bin/BayesIM-new # bin/imputeTwoTrait

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

bin/BayesIM:   $(OBJECTS_DIR)BayesIM.o	 $(OBJECTS_DIR)Configuration.o	  $(OBJECTS_DIR)utility.o	
	$(LINK) $(LFLAGS) -o bin/BayesIM  $(OBJECTS_DIR)BayesIM.o $(OBJECTS_DIR)Configuration.o $(OBJECTS_DIR)utility.o $(LIBS)

bin/BayesIMmt:   $(OBJECTS_DIR)BayesIMmt.o	 $(OBJECTS_DIR)Configuration.o	  $(OBJECTS_DIR)utility.o	
	$(LINK) $(LFLAGS) -o bin/BayesIMmt  $(OBJECTS_DIR)BayesIMmt.o $(OBJECTS_DIR)Configuration.o $(OBJECTS_DIR)utility.o $(LIBS)

bin/PullRegions:   $(OBJECTS_DIR)PullRegions.o	 $(OBJECTS_DIR)Configuration.o	  $(OBJECTS_DIR)utility.o	
	$(LINK) $(LFLAGS) -o bin/PullRegions  $(OBJECTS_DIR)PullRegions.o $(OBJECTS_DIR)Configuration.o $(OBJECTS_DIR)utility.o $(LIBS)




dist: 



clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) bin/BayesIM bin/PullRegions


 

####### Compile



$(OBJECTS_DIR)impute.o: $(SRC_DIR)impute.cpp $(HEADER_DIR)impute.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)impute.o $(SRC_DIR)impute.cpp

$(OBJECTS_DIR)imputeTwoTrait.o: $(SRC_DIR)imputeTwoTrait.cpp $(HEADER_DIR)impute.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)imputeTwoTrait.o $(SRC_DIR)imputeTwoTrait.cpp

$(OBJECTS_DIR)utility.o: $(SRC_DIR)utility.cpp $(HEADER_DIR)impute.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)utility.o $(SRC_DIR)utility.cpp 

$(OBJECTS_DIR)BayesIM.o: $(SRC_DIR)BayesIM.cpp $(HEADER_DIR)impute.h $(HEADER_DIR)Configuration.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)BayesIM.o $(SRC_DIR)BayesIM.cpp 

$(OBJECTS_DIR)BayesIMmt.o: $(SRC_DIR)BayesIMmt.cpp $(HEADER_DIR)impute.h $(HEADER_DIR)Configuration.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)BayesIMmt.o $(SRC_DIR)BayesIMmt.cpp 

$(OBJECTS_DIR)PullRegions.o: $(SRC_DIR)PullRegions.cpp $(HEADER_DIR)impute.h $(HEADER_DIR)Configuration.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)PullRegions.o $(SRC_DIR)PullRegions.cpp 

$(OBJECTS_DIR)BayesIM-new.o: $(SRC_DIR)BayesIM-new.cpp $(HEADER_DIR)impute.h $(HEADER_DIR)Configuration.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)BayesIM-new.o $(SRC_DIR)BayesIM-new.cpp 

$(OBJECTS_DIR)Configuration.o: $(SRC_DIR)Configuration.cpp  $(HEADER_DIR)Configuration.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)Configuration.o $(SRC_DIR)Configuration.cpp 

$(OBJECTS_DIR)imputeMCMC.o: $(SRC_DIR)imputeMCMC.cpp $(HEADER_DIR)impute.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)imputeMCMC.o $(SRC_DIR)imputeMCMC.cpp




####### Install

install:   $(BINARY_DIR)BayesIM   $(BINARY_DIR)BayesIMmt $(BINARY_DIR)PullRegions

$(BINARY_DIR)BayesIM: bin/BayesIM	
	cp bin/BayesIM $(BINARY_DIR)

$(BINARY_DIR)BayesIMmt: bin/BayesIMmt	
	cp bin/BayesIMmt $(BINARY_DIR)

$(BINARY_DIR)PullRegions: bin/PullRegions	
	cp bin/PullRegions $(BINARY_DIR) 

uninstall:   FORCE

FORCE:

