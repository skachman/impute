#############################################################################
# Makefile for building: GenSel
# on Darwin OS X
#############################################################################

####### Compiler, tools and options

MATVECLOC     = matvec
CC            = clang
CXX           = clang++
DEFINES       = 
CFLAGS        = -pipe -g -Wall -W $(DEFINES)
#CXXFLAGS      = -pipe -O3 -ftree-vectorize -ftree-vectorizer-verbose=5 -funsafe-math-optimizations -g  -arch x86_64 $(DEFINES) # -Wall -W
CXXFLAGS      = -pipe -O3 -funsafe-math-optimizations -g  -arch x86_64 $(DEFINES) #
#CXXFLAGS      = -pipe  -g  -arch x86_64 $(DEFINES) #-Wall -W

INCPATH       =  -I/opt/local/include/eigen3/ -I/opt/local/include -I. -Iinclude
LINK          = clang++
LFLAGS        = -prebind
LIBS          = $(SUBLIBS) /opt/local/lib/libmatvec.a 
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

####### Files

SOURCES       = $(SRC_DIR)impute.cpp  


OBJECTS       = $(OBJECTS_DIR)impute.o	

DESTDIR       = 

TARGET        = bin/impute 

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

all: Makefile $(TARGET) bin/imputeMCMC

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

bin/imputeMCMC:   $(OBJECTS_DIR)imputeMCMC.o	 
	$(LINK) $(LFLAGS) -o bin/imputeMCMC  $(OBJECTS_DIR)imputeMCMC.o $(OBJCOMP) $(LIBS)

dist: 
	@$(CHK_DIR_EXISTS) obj/GenSel1.0.0 || $(MKDIR) obj/GenSel1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) obj/GenSel1.0.0/ && $(COPY_FILE) --parents CMPBayesABC.cpp obj/GenSel1.0.0/ && (cd `dirname obj/GenSel1.0.0` && $(TAR) GenSel1.0.0.tar GenSel1.0.0 && $(COMPRESS) GenSel1.0.0.tar) && $(MOVE) `dirname obj/GenSel1.0.0`/GenSel1.0.0.tar.gz . && $(DEL_FILE) -r obj/GenSel1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_objective_c_make_all:
compiler_objective_c_clean:
compiler_moc_header_make_all:
compiler_moc_header_clean:
compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_rez_source_make_all:
compiler_rez_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile



$(OBJECTS_DIR)impute.o: $(SRC_DIR)impute.cpp $(HEADER_DIR)impute.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)impute.o $(SRC_DIR)impute.cpp

$(OBJECTS_DIR)imputeMCMC.o: $(SRC_DIR)imputeMCMC.cpp $(HEADER_DIR)impute.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(OBJECTS_DIR)imputeMCMC.o $(SRC_DIR)imputeMCMC.cpp




####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

