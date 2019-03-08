#include "Matrices.h"

#include "crism_mex_defines.h"

// functions to return pointers to Matrix?D objects
// WARNING: on heap, so dealloc them manually

Matrix1D wrapMatlabDoubleArray1D(const mxArray *matlabArray) {
    // just assumes proper dimensionality. All arrays will produce valid 1D output with this
    size_t num_els = mxGetNumberOfElements(matlabArray);
    
    Matrix1D matrix(num_els,false);
    matrix.setData((double*)mxGetData(matlabArray));
    matrix.setDataForfeit(false); // Don't destroy the input data when this matrix is destroyed
    
    return matrix;
}

Matrix2D wrapMatlabDoubleArray2D(const mxArray *matlabArray) {
    mwSize numDims = mxGetNumberOfDimensions(matlabArray);
    // Note: the docs for R2016a suggest that the size_t here should be mwSize, but in Matlab R2015a, it is size_t. Perhaps not a solid API.
    const size_t* dims = mxGetDimensions(matlabArray);
    if(numDims != 2) {
        printf("Improper number of dimensions in Array (%i instead of 2)\n",numDims);
        return Matrix2D(1,1);
    }
    
    size_t num_rows = dims[0];
    size_t num_cols = dims[1];
    
    Matrix2D matrix(num_rows,num_cols,false);
    matrix.setData((double*)mxGetData(matlabArray));
    matrix.setDataForfeit(false); // Don't destroy the input data when this matrix is destroyed
    
    return matrix;
}

Matrix3D wrapMatlabDoubleArray3D(const mxArray *matlabArray) {
    mwSize numDims = mxGetNumberOfDimensions(matlabArray);
    // Note: the docs for R2016a suggest that the size_t here should be mwSize, but in Matlab R2015a, it is size_t. Perhaps not a solid API.
    const size_t* dims = mxGetDimensions(matlabArray);
    if(numDims != 3) {
        printf("Improper number of dimensions in Array (%i instead of 3)\n",numDims);
        return Matrix3D(1,1,1);
    }
    
    size_t num_rows = dims[0];
    size_t num_cols = dims[1];
    size_t num_bands = dims[2];
    
    Matrix3D matrix(num_rows,num_cols,num_bands,false);
    matrix.setData((double*)mxGetData(matlabArray));
    matrix.setDataForfeit(false); // Don't destroy the input data when this matrix is destroyed
    
    return matrix;
}

double wrapMatlabDoubleScalar(const mxArray *matlabArray) {
    size_t numEls = mxGetNumberOfElements(matlabArray);
    if(numEls != 1) {
        printf("Matlab array does not represnt scalar! Passing back 0.\n");
        return 0.0;
    }
    
    double* data = (double*)mxGetData(matlabArray);
    return *data;
}