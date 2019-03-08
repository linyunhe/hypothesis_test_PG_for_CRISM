#ifndef __H_CRISM_MEX_DEFINES
#define __H_CRISM_MEX_DEFINES

#include "mex.h"

#ifdef MATLAB_MEX_FILE // Defined if this file is being compiled by MEX, which requires its own fancy print statement
    #define printf(...) mexPrintf( __VA_ARGS__ );mexEvalString("pause(.001);"); // the pause evaluation is necc. to flush the stream
#endif

#endif /* __H_CRISM_MEX_DEFINES */