#include <stdio.h>
#include "mex.h"

#ifdef _OPENMP
#include "omp.h"
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
	printf("Test mex function was successfully created!\n");
	printf("Number of inputs: %d\nNumber of outputs: %d\n",nrhs,nlhs);
	printf("Compiled on %s %s\n",__DATE__,__TIME__);
#ifdef __GNUC__
	printf("GCC version %d.%d.%d\n",__GNUC__, __GNUC_MINOR__,__GNUC_PATCHLEVEL__);
#endif
    
#ifdef _OPENMP
    printf("OpenMP is avaliable with %d cores\n", omp_get_num_procs());
#endif
    
}
