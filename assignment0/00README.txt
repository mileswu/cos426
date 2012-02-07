This directory contains skeleton code as a starting point for
Assignment 0 of COS 426. 


FILE STRUCTURE
==============

COMPILER.txt - notes on installing a C++ compiler and development environment.

Makefile - instructions for the "make" tool (used on Linux, Mac, and Cygwin)
 to run your program on the input images.

NMakefile - instructions for the "nmake" tool (used with Visual Studio)
 to run your program on the input images.

input/ - Contains example input images.

output/ - Empty to start.  It will contain images produced by your program.

output-sample/ - Examples of what your code should produce.

writeup.html - Skeleton HTML file for you to use as a basis for your writeup.

src/ - Directory with source code.
    Makefile - Unix/Mac/Cygwin makefile for building your code with "make".
    cos426_assignment0.sln - Project file for Visual Studio.
    imgpro.cpp - Main program: parses the command line arguments, and
     calls the appropriate image functions.
    R2Image.[cpp/h] - Image class with processing functions.
     ** This is the only file that you need to edit. **
    R2Pixel.[cpp/h] - Pixel class.
    morphlines.cpp - A program for specifying line correspondences between
     two images, which you will use for morphing.  You do not need to use
     it for Assignment 0, but make sure you can compile and execute it.
    R2/ - A library of useful 2D geometric primitives.  Built automatically.
    jpeg/ - A library for reading/writing JPEG files.  Built automatically.
    glut/ - The GLUT library for Visual Studio.


HOW TO PROCEED
==============

1. Read COMPILER.txt and install a C++ compiler on your development machine.

2. If you are using Visual Studio, open the .sln file, and build the solution.
Otherwise, cd into the "src" directory, and type "make".  STOP HERE AND ASK
FOR HELP if you cannot successfully compile the imgpro and morphlines programs.

3. Implement the R2Image::Brighten method in R2Image.cpp.  You will need to
read the R2Image.h and R2Pixel.h files to find out about the methods
implemented for these classes.  For this assignment (and only this
assignment), you may consult with fellow students or TAs to understand how
to implement the required function.

4. Recompile the code.

5. Run your code.  On Unix/Mac/Cygwin, change into the main assignment
directory (not the src/ subdirectory) and run "make".  With Visual Studio,
open the "Visual Studio Command Prompt", change into the main assignment
directory (not the src/ subdirectory) and run "nmake /f NMakefile".
Look for the output files in the output/ subdirectory.

6. Edit writeup.html to add a description of, and links to, the images
you just created.  Open writeup.html in a web browser to see the webpage.

7. Create a .zip file containing the contents of this directory using the 
naming convention that you can find on the Assignments page at
http://www.cs.princeton.edu/courses/archive/spr11/cos426/assignments.html#submitting

8. Submit it at the Dropbox link on the course webpage.

