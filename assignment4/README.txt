This directory contains skeleton code as a starting point for assignment4 of COS 426. 


FILE STRUCTURE
==============

There are several files, but you should mainly change src/particle.cpp.

  src/ - Directory with source code
    Makefile - makefile for building the project with "make".
    particleview.[vcproj/sln] - Project file for compiling particleview in Visual Studio on Windows
    particleview.cpp - Interactive program for viewing scenes
    particle.cpp - Main particle system simulation (this is the file you will edit)
    R3Scene.[cpp/h] - Class used for reading and storing scenes to be ray traced.
    R2/ - A library of useful 2D geometric primitives (includes R2Image and R2Pixel)
    R3/ - A library of useful 3D geometric primitives (includes R3Mesh)
    jpeg/ - A library for reading/writing JPEG files
  input/ - Contains example input scenes. 
  output/ - Is empty to start -- it will contain the movies and images produced by your program
  art/ - Empty to start -- it will (optionally) contain movies and art contest submissions produced by your program
  writeup.html - a skeleton HTML file that you can use as a basis for your writeup 


COMPILATION
===========

If you are developing on a Windows machine and have Visual Studio
installed, use the provided project solution files
(src/particleview.sln) to build the program. If you are developing on
a Mac, cd into the src/ directory and type "make".  If you are
developing on a Linux machine, cd into the src/ directory, then type
"cp Makefile.linux Makefile" and then you can type "make" to compile
the programs.  In any case, executables called particleview (or
particleview.exe) will be created in the src/ directory.
