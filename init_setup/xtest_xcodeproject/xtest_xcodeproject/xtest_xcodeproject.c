/************************************************************************************
Project xtest_xcodeproject created from Xcode template
modified by Fedor Iskhakov (fediskhakov@gmail.com, fedor.iskh.ru) 
http://www.mathworks.com.au/support/solutions/en/data/1-BOCDDE/index.html?solution=1-BOCDDE
************************************************************************************/
/*	
To debug using Xcode: 
 1. Compile in debug mode (with debug symbols)
 2. Attach to relevant Matlab process in [Product]->[Attach To Process]
 3. Place some breakpoints
 4. Run in Matlab as usual
 For OpenMP support:
 1. Install GCC on Xcode5 (current Apple LLVM 5.0 compiler does not support OpenMP)
 2. Add -fopenmp flag to compiler flags in mex.xcconfig file
 3. Add -fgomp in OtherLinkerFlags on the setting page of a particular target (upper left dropdown menu)
 */

#include "stdio.h"
#include "mex.h"

#ifdef _OPENMP
#include "omp.h"
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
	printf("Test mex function was successfully created!\n");
	printf("Number of inputs: %d\nNumber of outputs: %d\n",nrhs,nlhs);
	printf("Compiled on %s %s\n",__DATE__,__TIME__);
    
#ifdef _OPENMP
    printf("OpenMP is avaliable with %d cores\n", omp_get_num_procs());
#endif
    
}
