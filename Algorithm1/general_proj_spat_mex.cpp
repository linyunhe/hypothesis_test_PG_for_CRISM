#include "mex.h"
#include <omp.h>
#include <cmath>
#include "Matrices.h"
#include "MatricesMexWrappers.h"
#include "ompThreader.h"
#include "crism_mex_defines.h"

// Uses OpenMP. Compile from Matlab prompt with:
// mex -v general_proj_spat_mex.cpp -largeArrayDims CC=g++ COMPFLAGS="/openmp $COMPFLAGS" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" -Dcrism_omp_on


// constants via Matlab "format long"
#define PI 3.141592653589793
#define HALF_PI 1.570796326794897
#define PI_OVER_180 0.017453292519943

using namespace std;

void general_proj_spat_cpp(
	Matrix3D *in_matrix,     // FORWARD: c_num_rows x c_num_cols x num_bands; BACKWARD: num_rows x num_cols x num_bands
	bool isForward,
	bool isOffNadir,
	Matrix3D *ddr_iof,       // num_rows x num_cols x 14
	Matrix1D *ddr_lat,       // num_rows
	Matrix1D *ddr_lon,       // num_cols
	Matrix3D *sb,            // (1 or 10) x num_cols x num_bands
	int size_kernel,
	Matrix2D *lat_number,    // num_rows x num_cols					(contains indices into C: assumes they are 1-indexed)
	Matrix2D *lon_number,    // num_rows x num_cols					(contains indices into C: assumes they are 1-indexed)
	double pix_spacing,
	Matrix2D *mro_position,  // num_rows x 3
	Matrix3D *pix_positions, // c_num_rows x c_num_cols x 3
	Matrix3D *pix_normals,   // c_num_rows x c_num_cols x 3
	// output
	Matrix3D *conv_layer     // FORWARD: num_rows x num_cols x num_bands; BACKWARD: c_num_rows x c_num_cols x num_bands
)
{
	// NOTE: caller responsible for allocating the output matrix (and choosing its size)

	// NOTE: this version doesn't handle precomputation

	const int padding = (size_kernel - 1) / 2;
	const int num_elements = size_kernel * size_kernel;
	const double ifov = 61.5e-6;
	const double radius = 3396190.0;
	const double spectral_sampling = 6.55; //nm/channel
	const double nominal_altitude_sq = 266e3 * 266e3;

	// get dimensions
	size_t num_rows, num_cols, num_bands, c_num_rows, c_num_cols;
	if (isForward) {
		in_matrix->getDims(&c_num_rows, &c_num_cols, &num_bands);
		conv_layer->getDims(&num_rows, &num_cols, &num_bands);
	}
	else {
		in_matrix->getDims(&num_rows, &num_cols, &num_bands);
		conv_layer->getDims(&c_num_rows, &c_num_cols, &num_bands);
	}

	Matrix1D projected_areas(c_num_rows);


	if (isOffNadir) {
		double lon_spacing = ((*ddr_lon)(c_num_cols - 1) - (*ddr_lon)(0)) / (c_num_cols - 1);
		double lat_spacing = ((*ddr_lat)(c_num_rows - 1) - (*ddr_lat)(0)) / (c_num_rows - 1);
		for (int c_row = 0; c_row < c_num_rows; ++c_row) {
			double top_lat = (*ddr_lat)(c_row) - 0.5 * lat_spacing;
			double bottom_lat = (*ddr_lat)(c_row) + 0.5 * lat_spacing;
			projected_areas(c_row) = abs(lon_spacing*(cos(HALF_PI - top_lat) - cos(HALF_PI - bottom_lat)));
		}
	}
	else {
		projected_areas.fill(pix_spacing * pix_spacing);
	}

	// The main loop over all pixels
#pragma omp parallel for
	for (int col = 0; col < num_cols; ++col) {

		Matrix1D fwhm_alt1_sq(num_bands);
		Matrix1D leading_coeff(num_bands);

		Matrix1D LOS_target(3); // vector in 3-space (only needed for off-nadir case)
		double LOS_target_sq;

		for (int band = 0; band < num_bands; ++band) {
			fwhm_alt1_sq(band) = sb->ind(0, col, band) * ifov / spectral_sampling; // 0th row of sb contains spectral fwhms, which we now scale
			fwhm_alt1_sq(band) *= fwhm_alt1_sq(band);
		}

		for (int row = 0; row < num_rows; ++row) {
			// we know these arrays contain integral values
			int c_row = (int)lat_number->ind(row, col) - 1; // The "-1" is 1-based -> 0-based index conversion
			int c_col = (int)lon_number->ind(row, col) - 1; // The "-1" is 1-based -> 0-based index conversion

			double pix_lat, pix_lon;

			if (isOffNadir) {
				double R_target_i;

				LOS_target_sq = 0;
				for (int i = 0; i < 3; ++i) {
					R_target_i = pix_positions->ind(c_row, c_col, i);
					LOS_target.ind(i) = R_target_i - mro_position->ind(row, i);
					LOS_target_sq += LOS_target.ind(i) *  LOS_target.ind(i);

					// Omitting the filling of the center of apparent_areas that's here in the Matlab code, since the value is overwritten before being used
				}

			}
			else {
				LOS_target_sq = nominal_altitude_sq;
				pix_lat = ddr_iof->ind(row, col, 3)*PI_OVER_180; // Latitude slice (1-indexed, the 4th slice)
				pix_lon = ddr_iof->ind(row, col, 4)*PI_OVER_180; // Longitude slice (1-indexed, the 5th slice)
			}

			// build leading_coeff
			for (int band = 0; band < num_bands; ++band) {
				double fwhm_sq_b = fwhm_alt1_sq.ind(band) * LOS_target_sq;
				leading_coeff.ind(band) = log(16) / (PI * fwhm_sq_b);
			}

			Matrix3D fragment(size_kernel, size_kernel, num_bands); // not used if backprojection, but small enough to ignore waste
			Matrix2D apparent_areas(size_kernel, size_kernel);
			Matrix3D gaussian_kernel(size_kernel, size_kernel, num_bands);

			for (int row_inc = -padding; row_inc <= padding; ++row_inc) {

				int row_iter = c_row + row_inc;

				double iter_lat, y;

				// If current cell out of bounds, protected against access violations happening here,
				// but keep going anyway. It's cleaner to handle filling all the various arrays with the
				// default OOB value in only one place.
				if (row_iter >= 0 && row_iter < c_num_rows) {
					if (!isOffNadir) {
						iter_lat = ddr_lat->ind(row_iter);
						y = radius * (pix_lat - iter_lat);
					}
				}

				// will do this all at once inside
				//// out-of-bounds in rows
				//if(row_iter < 0 || row_iter > c_num_rows) {
				//    fragment = 
				//    continue;
				//}

				for (int col_inc = -padding; col_inc <= padding; ++col_inc) {

					int col_iter = c_col + col_inc;

					// out-of-bounds filling
					if (row_iter < 0 || row_iter >= c_num_rows || col_iter < 0 || col_iter >= c_num_cols) {
						apparent_areas.ind(row_inc + padding, col_inc + padding) = 1; // This will make the out-of-bounds values much lower
						for (int band = 0; band < num_bands; ++band) {
							if (isForward) {
								fragment.ind(row_inc + padding, col_inc + padding, band) = 0;
							}
							gaussian_kernel.ind(row_inc + padding, col_inc + padding, band) = leading_coeff.ind(band) * PI;
						}
						continue;
					}

					// within-bounds filling

					// for on-nadir
					double iter_lon, radius_eff, x;

					// for off-nadir
					Matrix1D R_neighbor(3), LOS_neighbor(3), normal_neighbor(3);
					double LOS_neighbor_sq = 0, cossq = 0, viewing_area_cosine = 0, local_radius_sq = 0, slope_cosine = 0;
					double tansq, ttn_dist_sq;

					if (isOffNadir) {
						for (int i = 0; i < 3; ++i) {
							R_neighbor.ind(i) = pix_positions->ind(row_iter, col_iter, i);
							LOS_neighbor.ind(i) = R_neighbor.ind(i) - mro_position->ind(row, i);
							LOS_neighbor_sq += LOS_neighbor.ind(i) * LOS_neighbor.ind(i);
							// find distance between target and projection of
							// neighbor onto plane intersecting target and
							// perpendicular to LOS_target
							cossq += LOS_target.ind(i)*LOS_neighbor.ind(i); // intermediate result

							// preliminaries to populate apparent area matrix
							normal_neighbor.ind(i) = pix_normals->ind(row_iter, col_iter, i);
							viewing_area_cosine += LOS_neighbor.ind(i)*normal_neighbor.ind(i);
							local_radius_sq += R_neighbor.ind(i)*R_neighbor.ind(i);
							slope_cosine += R_neighbor.ind(i)*normal_neighbor.ind(i);
						}

						cossq = (cossq * cossq) / (LOS_target_sq * LOS_neighbor_sq); // calculate final value (numerator is sum from last loop)
						tansq = (1 / cossq) - 1;
						ttn_dist_sq = LOS_target_sq*tansq;

						viewing_area_cosine = abs(viewing_area_cosine / sqrt(LOS_neighbor_sq)); // calculate final value (one term is sum from last loop)
						slope_cosine = abs(slope_cosine / sqrt(local_radius_sq)); // calculate final value (one term is sum from last loop)

						// populate apparent area matrix
						apparent_areas.ind(row_inc + padding, col_inc + padding) = local_radius_sq * projected_areas(row_iter) * viewing_area_cosine / slope_cosine;

					}
					else {
						iter_lon = ddr_lon->ind(col_iter);
						radius_eff = radius*sin(HALF_PI - abs(pix_lat + iter_lat) / 2.0);
						x = radius_eff*(pix_lon - iter_lon);
						apparent_areas.ind(row_inc + padding, col_inc + padding) = projected_areas(row_iter);
					}

					for (int band = 0; band < num_bands; ++band) {
						double leading_coeff_i = leading_coeff.ind(band);
						if (isOffNadir) {
							// compute the (normalized) gaussian:
							gaussian_kernel.ind(row_inc + padding, col_inc + padding, band) = leading_coeff_i * exp(-ttn_dist_sq*PI*leading_coeff_i);
						}
						else {
							gaussian_kernel.ind(row_inc + padding, col_inc + padding, band) = leading_coeff_i * exp(-PI*leading_coeff_i*(x*x + y*y));
						}

						if (isForward) {
							fragment.ind(row_inc + padding, col_inc + padding, band) = in_matrix->ind(row_iter, col_iter, band);
						}
					}
				}
			}

			// Multiply Gaussian kernel by the apparent area matrix
			for (int c = 0; c < size_kernel; ++c) {
				for (int r = 0; r < size_kernel; ++r) {
					for (int b = 0; b < num_bands; ++b) {
						gaussian_kernel.ind(r, c, b) = gaussian_kernel.ind(r, c, b) * apparent_areas.ind(r, c);
					}
				}
			}

			// Normalizing the kernel (in off-nadir case)
			if (isOffNadir) {
				for (int b = 0; b < num_bands; ++b) {
					double gaussian_kernel_volume_i = 0;
					for (int c = 0; c < size_kernel; ++c) {
						for (int r = 0; r < size_kernel; ++r) {
							gaussian_kernel_volume_i += gaussian_kernel.ind(r, c, b);
						}
					}
					// now volume is in hand; divide it through all elements in the band
					for (int c = 0; c < size_kernel; ++c) {
						for (int r = 0; r < size_kernel; ++r) {
							gaussian_kernel.ind(r, c, b) /= gaussian_kernel_volume_i;
						}
					}
				}
			}

			// Convolution
			// Forward Projection Case: Dot product of gaussian_kernel and fragment per band, stored for current row & col
			if (isForward) {
				for (int b = 0; b < num_bands; ++b) {
					double dot_result = 0;
					for (int c = 0; c < size_kernel; ++c) {
						for (int r = 0; r < size_kernel; ++r) {
							dot_result += gaussian_kernel.ind(r, c, b) * fragment.ind(r, c, b);
						}
					}
					conv_layer->ind(row, col, b) = dot_result;
				}
			}
			else {

				for (int band = 0; band < num_bands; ++band) {
					for (int c = 0; c < size_kernel; ++c) { // column in iterating over kernel
						int col_in_scene = c + c_col - padding;

						for (int r = 0; r < size_kernel; ++r) { // row in iterating over kernel
							int row_in_scene = r + c_row - padding;

							// out-of-bounds check: if so, skip.
							if (row_in_scene < 0 || row_in_scene >= c_num_rows || col_in_scene < 0 || col_in_scene >= c_num_cols) {
								continue;
							}

							double output_increment = in_matrix->ind(row, col, band) * gaussian_kernel.ind(r, c, band);

#pragma omp atomic
							conv_layer->ind(row_in_scene, col_in_scene, band) += output_increment;

						}
					}

				}

			}



		}

	}

	if (isForward) {
		printf("Spatial Forward Projection (C++):\n");
	}
	else {
		printf("Spatial Backprojection (C++):\n");
	}

}

/* 
 * A MEX version of calculate_spatial_kernel_anygeo.m, to replace its forward
 * and backward projection modes (except for the version that stores a
 * precomputed kernel)
 *
 
 
 * Inputs:
 * 0. in_matrix
 * 1. isOffNadir
 * 2. ddr_iof
 * 3. ddr_lat
 * 4. ddr_lon
 * 5. sb
 * 6. size_kernel
 * 7. lat_number
 * 8. lon_number
 * 9. pix_spacing
 * 10. mro_position
 * 11. pix_positions
 * 12. pix_normals
 
 * Outputs:
 * 0. conv_layer
 *
 */

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    #define in_matrix_in      prhs[0]
	#define isForward_in      prhs[1]
    #define isOffNadir_in     prhs[2]
    #define ddr_iof_in        prhs[3]
    #define ddr_lat_in        prhs[4]
    #define ddr_lon_in        prhs[5]
    #define sb_in             prhs[6]
    #define size_kernel_in    prhs[7]
    #define lat_number_in     prhs[8]
    #define lon_number_in     prhs[9]
    #define pix_spacing_in    prhs[10]
    #define mro_position_in   prhs[11]
    #define pix_positions_in  prhs[12]
    #define pix_normals_in    prhs[13]
	#define omp_threads_in    prhs[14]
    
    #define conv_layer_out   plhs[0]
    
    // Checking the properties of the inputs & outputs
    // TODO?
    
    // Wrap the input Matlab matrices
    // (Allocates on stack inside and uses default copy constructor to get it out here)
    Matrix3D in_matrix = wrapMatlabDoubleArray3D(in_matrix_in);
	bool isForward = wrapMatlabDoubleScalar(isForward_in) != 0.0; // scalar double -> bool conversion
    bool isOffNadir = wrapMatlabDoubleScalar(isOffNadir_in) != 0.0; // scalar double -> bool conversion
    Matrix3D ddr_iof = wrapMatlabDoubleArray3D(ddr_iof_in);
    Matrix1D ddr_lat = wrapMatlabDoubleArray1D(ddr_lat_in);
    Matrix1D ddr_lon = wrapMatlabDoubleArray1D(ddr_lon_in);
    Matrix3D sb = wrapMatlabDoubleArray3D(sb_in);
    int size_kernel = (int)wrapMatlabDoubleScalar(size_kernel_in);
    Matrix2D lat_number = wrapMatlabDoubleArray2D(lat_number_in);
    Matrix2D lon_number = wrapMatlabDoubleArray2D(lon_number_in);
    double pix_spacing = wrapMatlabDoubleScalar(pix_spacing_in);
    // only give serious values to these three if we're offnadir (otherwise, dimension check fails)
    Matrix2D mro_position = isOffNadir ? wrapMatlabDoubleArray2D(mro_position_in) : Matrix2D(1,1,false);
    Matrix3D pix_positions = isOffNadir ? wrapMatlabDoubleArray3D(pix_positions_in) : Matrix3D(1,1,1,false);
    Matrix3D pix_normals = isOffNadir ? wrapMatlabDoubleArray3D(pix_normals_in) : Matrix3D(1,1,1,false);
    
    // output matrix: determining its dimensions
    size_t num_rows, num_cols, num_bands, c_num_rows, c_num_cols;
    in_matrix.getDims(nullptr,nullptr,&num_bands);
	ddr_lat.getDims(&c_num_rows);
	ddr_lon.getDims(&c_num_cols);
    lat_number.getDims(&num_rows,&num_cols);
    
	size_t output_dims[3];
	if (isForward) {
		output_dims[0] = num_rows;
		output_dims[1] = num_cols;
		output_dims[2] = num_bands;
	}
	else {
		output_dims[0] = c_num_rows;
		output_dims[1] = c_num_cols;
		output_dims[2] = num_bands;
	}
    
    // setup of output Matrix somewhere Matlab can get to it
	conv_layer_out = mxCreateNumericArray(3, output_dims, mxDOUBLE_CLASS, mxREAL);
    Matrix3D conv_layer(output_dims[0], output_dims[1], output_dims[2], false);
    conv_layer.setData(mxGetPr(conv_layer_out)); // use Matlab-alloc'd space to let us carry result out
    conv_layer.setDataForfeit(false); // allows us to carry out our data

	 // setup OpenMP
	int omp_threads = std::floor(mxGetScalar(omp_threads_in));
	setupOpenMP(omp_threads);
    
    // Call the c++ function
	general_proj_spat_cpp(&in_matrix,isForward,isOffNadir,&ddr_iof,&ddr_lat,&ddr_lon,&sb,size_kernel,&lat_number,&lon_number,pix_spacing,&mro_position,&pix_positions,&pix_normals,&conv_layer);
    
}

