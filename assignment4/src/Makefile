#
# Mac makefile for assignment #4
#

# 
# List of source files
#

SRCS=particleview.cpp particle.cpp R3Scene.cpp raytrace.cpp
OBJS=$(SRCS:.cpp=.o)

#
# Compile and link options
#

CXX=g++
CXXFLAGS=-Wall -I. -g -DUSE_JPEG -O2

#
# Libraries
#

LIBS=R3/libR3.a R2/libR2.a jpeg/libjpeg.a

#
# OpenGL libraries
#
UNAME := $(shell uname)
ifneq (,$(findstring Darwin,$(UNAME)))
	GLLIBS = -framework GLUT -framework OpenGL
else
  ifneq (,$(findstring CYGWIN,$(UNAME)))
	GLLIBS = -lglut32 -lglut -lGLU -lGL
  else
	GLLIBS = -lglut -lGLU -lGL
  endif
endif

#
# GNU Make: targets that don't build files
#

.PHONY: all clean distclean


#
# Make targets
#

all: particleview

particleview: $(OBJS) $(LIBS)
		rm -f $@
	    $(CXX) $(CXXFLAGS) $^ -lm -o $@ $(GLLIBS)

R3/libR3.a: 
	    $(MAKE) -C R3

R2/libR2.a: 
	    $(MAKE) -C R2

jpeg/libjpeg.a: 
	    $(MAKE) -C jpeg

clean:
	    rm -f *.o particleview
		$(MAKE) -C R3 clean
		$(MAKE) -C R2 clean
		$(MAKE) -C jpeg clean

distclean:  clean
	    ${RM} -f *~ 

