/*
 * Matrices
 *
 * based on class found at:
 * http://www.cplusplus.com/forum/articles/17108/
 *
 * AN IMPORTANT DETAIL different from usual C matrices:
 * The elements are stored internally in the order that MATLAB uses.
 * This is only visible when using these classes if:
 *  1) You are importing/exporting flat arrays into these matrices
 *  2) You are indexing into them 1-dimensionally.
 *
 */

#include <iostream>
#include <cassert>

#include "crism_mex_defines.h"

#ifndef __H_MATRICES
#define __H_MATRICES

using namespace std;

class Matrix1D {
private:

	size_t nRows;
	double *arr;
    bool isDataForfeit;
    
public:

	// constructor
	Matrix1D(size_t rows, bool allocateSpace = true)
		: nRows(rows), arr(0), isDataForfeit(true) {

        if(allocateSpace) {
            initializeData();
        }
	}

	// destructor
	~Matrix1D()
	{
        if(isDataForfeit) {
            delete[] arr;
        }
	}

	// indexing
	double& ind(size_t r) {
        assert(r >= 0 && r < nRows);
		return arr[r];
	}

	// simpler 1-D indexing too
	double& operator() (size_t i) {
        assert(i >= 0 && i < nRows);
		return arr[i];
	}

	// get dims
	size_t num_rows() { return nRows; }
	size_t num_els() { return nRows; }
	void getDims(size_t *rows) {
		if (rows) {
			*rows = nRows;
		}
	}

	// ways to reset entire matrix to a value

	void fill(double val) {
		for (size_t i = 0; i < nRows; ++i) {
			arr[i] = val;
		}
	}

	void zeros() { fill(0); }
	void ones() { fill(1); }

	// print matrix to display
	void print(ostream& outStream) {
		// test getting it with ind()
		for (size_t i = 0; i < nRows; ++i) {
			outStream << this->ind(i) << "\t";
		}
        outStream << endl;
	}
    
    // print matrix with printf
	void print() {
		// test getting it with ind()
		for (size_t i = 0; i < nRows; ++i) {
			printf("%f\t",this->ind(i));
		}
        printf("\n");
	}
    
        
    // gets the pointer to the internal data. Does not make a copy!
    double* getData() {
        return arr;
    }
    
    // sets the pointer to the internal data. Does not make a copy!
    void setData(double *data) {
        arr = data;
    }
    
    // allocates space on heap for this matrix
    void initializeData() {
        if (nRows > 0) {
			arr = new double[nRows];
		}
    }
    
    // choose whether the internal data of this array will be deleted when it is destructed.
    // The default value when this has not been called is true.
    void setDataForfeit(bool isDataForfeit) {
        this->isDataForfeit = isDataForfeit;
    }

};

class Matrix2D {
private:
	
	size_t nRows, nCols;
	double *arr;
    bool isDataForfeit;
    
public:

	// constructor
	Matrix2D(size_t rows, size_t cols, bool allocateSpace = true)
		: nRows(rows), nCols(cols), arr(0), isDataForfeit(true) {

        if(allocateSpace) {
            initializeData();
        }
	}

	// destructor
	~Matrix2D()
	{
        if(isDataForfeit) {
            delete[] arr;
        }
	}

	// indexing
	double& ind(size_t r, size_t c) {
        assert(r >= 0 && r < nRows && c >= 0 && c < nCols);
		return arr[c*nRows + r];
	}

	// 1-D indexing. For ease of use, since we're using it so much
	double& operator() (size_t i) {
		return arr[i];
	}

	// get dims
	size_t num_rows() { return nRows; }
	size_t num_cols() { return nCols; }
	size_t num_els() { return nRows * nCols; }
	void getDims(size_t *rows, size_t *cols) {
		if (rows) {
			*rows = nRows;
		}

		if (cols) {
			*cols = nCols;
		}
	}

	// ways to reset entire matrix to a value

	void fill(double val) {
		for (size_t i = 0; i < nRows*nCols; ++i) {
			arr[i] = val;
		}
	}

	void zeros() { fill(0); }
	void ones() { fill(1); }

	// print matrix to display
	void print(ostream& outStream) {
		// test getting it with ind()
		for (size_t i = 0; i < nRows; ++i) {
			for (size_t j = 0; j < nCols; ++j) {
				outStream << this->ind(i, j) << "\t";
			}
			outStream << endl;
		}
	}
    
    // print matrix with printf
	void print() {
		// test getting it with ind()
		for (size_t i = 0; i < nRows; ++i) {
			for (size_t j = 0; j < nCols; ++j) {
				printf("%f\t",this->ind(i, j));
			}
			printf("\n");
		}
	}
        
    // gets the pointer to the internal data. Does not make a copy!
    double* getData() {
        return arr;
    }
    
    // sets the pointer to the internal data. Does not make a copy!
    void setData(double *data) {
        arr = data;
    }
    
    // allocates space on heap for this matrix
    void initializeData() {
        if (nRows > 0 && nCols > 0) {
			arr = new double[nRows * nCols];
		}
    }
    
    // choose whether the internal data of this array will be deleted when it is destructed.
    // The default value when this has not been called is true.
    void setDataForfeit(bool isDataForfeit) {
        this->isDataForfeit = isDataForfeit;
    }

};

class Matrix3D {
private:

	size_t nRows, nCols, nBands;
	double *arr;
    bool isDataForfeit;
    
public:

	// constructor  
    Matrix3D(size_t rows, size_t cols, size_t bands, bool allocateSpace = true)
		: nRows(rows), nCols(cols), nBands(bands), arr(0), isDataForfeit(true) {
        
        if(allocateSpace) {
            initializeData();
        }
	}

	// destructor
	~Matrix3D()
	{
        if(isDataForfeit) {
            delete[] arr;
        }
	}

	// indexing
	double& ind(size_t r, size_t c, size_t b) {
        assert(r >= 0 && r < nRows && c >= 0 && c < nCols && b >= 0 && b < nBands);
		return arr[b*nRows*nCols + c*nRows + r];
	}
	// alias of 3d indexing
	double& operator() (size_t r, size_t c, size_t b) {
        assert(r >= 0 && r < nRows && c >= 0 && c < nCols && b >= 0 && b < nBands);
		return arr[b*nRows*nCols + c*nRows + r];
	}
	// 1-D indexing. For ease of use, since we're using it so much
	double& operator() (size_t i) {
		return arr[i];
	}

	// get dims
	size_t num_rows() { return nRows; }
	size_t num_cols() { return nCols; }
	size_t num_bands() { return nBands; }
	size_t num_els() { return nRows * nCols * nBands; }
	void getDims(size_t *rows, size_t *cols, size_t *bands) {
		if (rows) {
			*rows = nRows;
		}
		if (cols) {
			*cols = nCols;
		}
		if (bands) {
			*bands = nBands;
		}
	}

	// ways to reset entire matrix to a value

	void fill(double val) {
		for (size_t i = 0; i < nRows*nCols*nBands; ++i) {
			arr[i] = val;
		}
	}

	void zeros() { fill(0); }
	void ones() { fill(1); }

	// print matrix to display
	void print(ostream& outStream) {
		for (size_t k = 0; k < nBands; ++k) {
			for (size_t i = 0; i < nRows; ++i) {
				for (size_t j = 0; j < nCols; ++j) {
					outStream << this->ind(i, j, k) << "\t";
				}
				outStream << endl;
			}
			outStream << endl; // extra newline to separate slices in display
		}
	}

    // print matrix with printf
    void print() {
		for (size_t k = 0; k < nBands; ++k) {
			for (size_t i = 0; i < nRows; ++i) {
				for (size_t j = 0; j < nCols; ++j) {
					printf("%f\t",this->ind(i, j, k));
				}
                printf("\n");
			}
            printf("\n"); // extra newline to separate slices in display
		}
	}
    
    // gets the pointer to the internal data. Does not make a copy!
    double* getData() {
        return arr;
    }
    
    // sets the pointer to the internal data. Does not make a copy!
    void setData(double *data) {
        arr = data;
    }
    
    // allocates space on heap for this matrix
    void initializeData() {
        if (nRows > 0 && nCols > 0 && nBands > 0) {
			arr = new double[nRows * nCols * nBands];
		}
    }
    
    // choose whether the internal data of this array will be deleted when it is deconstructed.
    // The default value when this has not been called is true.
    void setDataForfeit(bool isDataForfeit) {
        this->isDataForfeit = isDataForfeit;
    }

};

class Matrix3D_bool {
	// TODO: might it be more efficient to manually pack the bits into ints or something?
    // WARNING: UNLIKE THE DOUBLE MATRICES ABOVE, THIS ONE HAS DIFFERENT SEMANTICS RELATED TO DELETING ITS CONTENTS UPON DECONSTRUCTION
private:

	size_t nRows, nCols, nBands;
	bool *arr;
public:

	// constructor
	Matrix3D_bool(size_t rows, size_t cols, size_t bands)
		: nRows(rows), nCols(cols), nBands(bands), arr(0) {

		if (nRows > 0 && nCols > 0 && nBands > 0) {
			arr = new bool[nRows * nCols * nBands];
		}
	}

	// destructor
	~Matrix3D_bool()
	{
		delete[] arr;
	}

	// indexing
	bool& ind(size_t r, size_t c, size_t b) {
		return arr[b*nRows*nCols + c*nRows + r];
	}
	// 1-D indexing. For ease of use, since we're using it so much
	bool& operator() (size_t i) {
		return arr[i];
	}

	// get dims
	size_t num_rows() { return nRows; }
	size_t num_cols() { return nCols; }
	size_t num_bands() { return nBands; }
	size_t num_els() { return nRows * nCols * nBands; }
	void getDims(size_t *rows, size_t *cols, size_t *bands) {
		if (rows) {
			*rows = nRows;
		}
		if (cols) {
			*cols = nCols;
		}
		if (bands) {
			*bands = nBands;
		}
	}

	// ways to reset entire matrix to a value

	void fill(bool val) {
		for (size_t i = 0; i < nRows*nCols*nBands; ++i) {
			arr[i] = val;
		}
	}

	void falses() { fill(false); }
	void trues() { fill(true); }

	// print matrix to display
	void print(ostream& outStream) {
		for (size_t k = 0; k < nBands; ++k) {
			for (size_t i = 0; i < nRows; ++i) {
				for (size_t j = 0; j < nCols; ++j) {
					outStream << this->ind(i, j, k) << "\t";
				}
				outStream << endl;
			}
			outStream << endl; // extra newline to separate slices in display
		}
	}
    
    // print matrix with printf
    void print() {
		for (size_t k = 0; k < nBands; ++k) {
			for (size_t i = 0; i < nRows; ++i) {
				for (size_t j = 0; j < nCols; ++j) {
					printf("%f\t",this->ind(i, j, k));
				}
                printf("\n");
			}
            printf("\n"); // extra newline to separate slices in display
		}
	}

};

#endif /* __H_MATRICES */
