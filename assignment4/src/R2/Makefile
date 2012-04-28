#
# Makefile for R2 library.
#

CXX=g++
CXXFLAGS=-Wall -I. -I.. -g -DUSE_JPEG

ifneq (,$(findstring Darwin,$(shell uname)))
  AR=libtool -static -o
else
  AR=ar -rcs
endif

#
# Output library name
#

LIB=./libR2.a

#
# List of files
#

CCSRCS=R2Distance.cpp \
  R2Segment.cpp R2Line.cpp R2Vector.cpp R2Point.cpp \
  R2Image.cpp R2Pixel.cpp
OBJS=$(CCSRCS:.cpp=.o)
INCS=$(HSRCS) $(CCSRCS:.cpp=.h)

$(LIB): $(OBJS)
	    - rm -f $(LIB)
	    $(AR) $@ $^

clean:
	    -  rm -f *~ *.o *.obj *.pdb $(LIB)


