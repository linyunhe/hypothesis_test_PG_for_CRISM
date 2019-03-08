#include "mex.h"

// A very simple check to see if our MEX C++ binaries can be run on this system or not

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    plhs[0] = mxCreateDoubleScalar(1);
}
