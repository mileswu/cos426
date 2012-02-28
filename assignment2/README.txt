This directory contains skeleton code as a starting point for assignment2 of COS 426. 


FILE STRUCTURE
==============

There are several files, but you should mainly change src/R3Mesh.cpp.

  src/ - Directory with source code
    Makefile - Unix/Mac makefile for building the project with "make". 
    meshpro.[vcproj/sln/suo] - Project file for Visual Studio 2005 on Windows
    meshpro.cpp - Main program, parses the command line arguments, and calls the appropriate mesh processing functions
    meshview.cpp - Interactive program for viewing meshes and making images of meshes
    R3Mesh.[cpp/h] - Mesh class with processing functions (this is the only file that you need to edit)
    R2/ - A library of useful 2D geometric primitives (includes R2Image and R2Pixel)
    R3/ - A library of useful 3D geometric primitives
    jpeg/ - A library for reading/writing JPEG files
  input/ - Contains example input images. 
  output/ - Empty to start -- it will contain the images produced by your program (see below)
  art/ - Empty to start -- it will (optionally) contain movies and art contest submissions produced by your program
  runme.bat - a script (for Windows) that you will fill in to demonstrate execution of your program
  runme.sh - same as <code>runme.bat, but for Mac OS X
  writeup.html - a skeleton HTML file that you can use as a basis for your writeup 


COMPILATION
===========

If you are developing on a Windows machine and have Visual Studio
installed, use the provided project solution file (assn2.sln) in the
src/ directory to build the program. If you are developing on a Mac or
Linux machine, cd into the src/ directory and type "make". In either
case, executables called meshpro and meshview (or meshpro.exe and meshview.exe) 
will be created in the src/ directory.
