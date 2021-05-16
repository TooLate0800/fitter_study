
MAKEFILE      = Makefile
####### Compiler, tools and options

CC            = gcc
CXX           = g++
CXXFLAGS      = -pipe -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -Wno-missing-field-initializers -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -O2 -std=gnu++0x -Wall -W -D_REENTRANT -fPIC $(DEFINES)
INCPATH       = -I. -I${ROOTSYS}/include -I$(PRAD_PATH)/include
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
LINK          = g++
LFLAGS        = -Wl,-O1 -Wl,-z,relro
LIBS          = -L${ROOTSYS}/lib -lCore -lRint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui -lSpectrum -lpthread -Llib 
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Files
              
TARGET        = RobustGenerator

first: all

####### Build rules
all: Makefile $(TARGET)

clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) $(TARGET)


RobustGenerator: RobustGenerator_loop_cut.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCPATH) $(LIBS)

install:   FORCE

uninstall:   FORCE

FORCE:


