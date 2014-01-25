EGDST SOFTWARE

ZICE2014 release

By: Fedor Iskhakov, fediskhakov@gmail.com http://fedor.iskh.ru

______________________________________________________________________________

This collection of files contains the latest version of the EGDST software and
some code snippets prepared for  the ZICE 2014.

In particular:

1. Directory "init_setup" contains some instructions and simple examples that
help to establish if C compiler and     other components are properly
installed in the system so that it can be used to run EGDST 2. Directory
"matlab_examples" contain some Matlab code to illustrate the principles of
Endogenous Gridpoint Methods 3. Directory "@egdstmodel" contains all the files
of the EGDST software package 4. Directory "egdst_examples" contains the
examples of use of EGDST software


_________________________ INSTALLATION _______________________________________

EGDST software requires C files to be compiled on the system, therefore the
first step is to install C compiler and make sure it works. Installation of C
compiler is specific to the system, see
http://www.mathworks.com.au/help/matlab/ref/mex.html

Folder "init_setup" contains a simple C code file test.c for testing the mex
setup.  Browse into this folder and run mex test.c test OR simply start

To enable OpenMP support, make sure you are using OpenMP supporting compiler
and to test the environment run mex text.c CFLAGS='\$CFLAGS -fopenmp'
LDFLAGS='\$LDFLAGS -fopenmp' test OR simply start

The output should mension how many cores are available for OpenMP, if the
output is identical to the run above  OpenMP is not enabled.


______________________________ RUNNING _______________________________________

Browse into "egdst_examples" folder, run "start" from Matlab

