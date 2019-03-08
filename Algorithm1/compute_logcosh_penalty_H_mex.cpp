#include "mex.h"
#include <omp.h>
#include <iostream>
#include "Matrices.h"
#include "PenaltyParams.h"
#include "ompThreader.h"
#include "crism_mex_defines.h"

// Uses OpenMP. Compile from Matlab prompt with:
// mex -v compute_new_image_logcosh_prior_mex.cpp -largeArrayDims CC=g++ COMPFLAGS="/openmp $COMPFLAGS" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" -Dcrism_omp_on

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);

/*
 * ---------------------------------------------------------------------------------------------------------------
 * -- COMPUTE_ORIGINAL_PENALTY_FUNCTION_SPAT_CPP: the C++ equivalent to compute_original_penalty_function_spat. --
 * ---------------------------------------------------------------------------------------------------------------
 */
void compute_original_penalty_function_spat_cpp(
	Matrix3D *sensitivity_image,
	Matrix3D *new_estimate,
	PenaltyParams penalty_params,
	// outputs
	Matrix3D *penalty_function
	)
{
    penalty_function->zeros(); // initialize
    
    size_t num_rows, num_cols, num_bands;
	sensitivity_image->getDims(&num_rows, &num_cols, &num_bands);
    
    // Spatial transforms
    if(penalty_params.spat_penalty_flag) {
        for(int row_inc = -1; row_inc <= 1; ++row_inc) {
            for(int col_inc = -1; col_inc <= 1; ++col_inc) {

                // skip the one iter of this we don't care about
                if ( row_inc == 0 && col_inc == 0 ) {
                    continue;
                }

                double neighbor_weight = 1 / sqrt( (double)row_inc * row_inc + col_inc * col_inc);

                // for each spatial point, we consider (r,c,:) of new_est as the curr. pixel and (r-row_inc,c-col_inc,:) of image_last_iter as the neighbor
                #pragma omp parallel for
                for (int r = 0; r < num_rows; ++r) {
                    int r_neigh = r + row_inc;
                    if (r_neigh < 0 || r_neigh >= num_rows) {
                        continue; // neighbor is in row w/o correspondance to any row in image
                    }

                    for (int c = 0; c < num_cols; ++c) {

                        int c_neigh = c + col_inc;
                        if (c_neigh < 0 || c_neigh >= num_cols) {
                            continue; // neighbor is in col w/o correspondance to any col in image
                        }
                        // now we can stop worrying about this point doing any improper indexing

                        // If neighbor_mask would be false, bail on this spatial point

                        double sensitivity_ij = (*sensitivity_image)(r, c, 0); // for spatial work, all bands at a point will be valid if any are
                        double sensitivity_neigh_ij = (*sensitivity_image)(r_neigh, c_neigh, 0); // for spatial work, all bands at a point will be valid if any are

                        // make sure we don't increment penalty_* for any point for which neighbour_mask would be false
                        if (sensitivity_ij == 0 || sensitivity_ij != sensitivity_ij) { // if this cell zero or NaN, don't bother
                            continue;
                        }

                        if (sensitivity_neigh_ij == 0 || sensitivity_neigh_ij != sensitivity_neigh_ij) { // if neighbor zero or NaN, don't bother
                            continue;
                        }

                        for (int b = 0; b < num_bands; ++b) {
                            // Increment the pixel corresonding to this (r,c,b) point's penalty_* values as appropriate

                            // bracketterm = 2*new_estimate_padded - c_last_iter_padded - temp;
                            double bracketterm_ijk = (*new_estimate)(r, c, b) - (*new_estimate)(r_neigh, c_neigh, b);
                            double trigArg_ijk = bracketterm_ijk / penalty_params.delta_spat;

                            double cosh_ijk = cosh(trigArg_ijk);

                            // updating penalty_function
                            if (trigArg_ijk > 500) { // approximate log-cosh by abs value
                                (*penalty_function)(r, c, b) += neighbor_weight * penalty_params.beta_spat * penalty_params.delta_spat * penalty_params.delta_spat * trigArg_ijk;
                            }
                            else { // use actual log-cosh
                                (*penalty_function)(r, c, b) += neighbor_weight * penalty_params.beta_spat * penalty_params.delta_spat * penalty_params.delta_spat * log(cosh_ijk);
                            }
                        }
                    }
                }

            }
        }
    }

}

/*
* ---------------------------------------------------------------------------------------------------------------
* -- COMPUTE_ORIGINAL_PENALTY_FUNCTION_SPEC_CPP: the C++ equivalent to compute_original_penalty_function_spec. --
* ---------------------------------------------------------------------------------------------------------------
*/
void compute_original_penalty_function_spec_cpp(
	Matrix3D *sensitivity_image,
	Matrix3D *new_estimate,
	PenaltyParams penalty_params,
	// outputs
	Matrix3D *penalty_function
)
{
	penalty_function->zeros(); // initialize

	size_t num_rows, num_cols, num_bands;
	sensitivity_image->getDims(&num_rows, &num_cols, &num_bands);

	// Spectral transforms
	if (penalty_params.spec_penalty_flag) {
		int span = penalty_params.spectral_span;
		for (int band_inc = -span; band_inc <= span; ++band_inc) {
			// skip the null transformation
			if (band_inc == 0) {
				continue;
			}

			double neighbor_weight = 1 / (double)abs(band_inc);

#pragma omp parallel for
			for (int b = 0; b < num_bands; ++b) {

				int b_neigh = b + band_inc;
				if (b_neigh < 0 || b_neigh >= num_bands) {
					continue; // neighbor is in row w/o correspondance to any row in image
				}

				for (int r = 0; r < num_rows; ++r) {
					for (int c = 0; c < num_cols; ++c) {

						// If neighbor_mask would be false, bail on this spatial point

						double sensitivity_ij = (*sensitivity_image)(r, c, 0);

						// make sure we don't increment penalty_* for any spatial location where the sensitivity image would be zero
						if (sensitivity_ij == 0 || sensitivity_ij != sensitivity_ij) {
							continue;
						}

						// Increment the pixel corresonding to this (r,c,b) point's penalty_* values as appropriate

						// bracketterm = 2*new_estimate_padded - c_last_iter_padded - temp;
						double bracketterm_ijk = (*new_estimate)(r, c, b) - (*new_estimate)(r, c, b_neigh);
						double trigArg_ijk = bracketterm_ijk / penalty_params.delta_spec;

						double cosh_ijk = cosh(trigArg_ijk);

						// updating penalty_function
						if (trigArg_ijk > 500) { // approximate log-cosh by abs value
							(*penalty_function)(r, c, b) += neighbor_weight * penalty_params.beta_spec * penalty_params.delta_spec * penalty_params.delta_spec * trigArg_ijk;
						}
						else { // use actual log-cosh
							(*penalty_function)(r, c, b) += neighbor_weight * penalty_params.beta_spec * penalty_params.delta_spec * penalty_params.delta_spec * log(cosh_ijk);
						}
					}
				}
			}
		}
	}
}

/*
 * ---------------------------------------------------------------------------------------------------------
 * -- COMPUTE_DECOMPOSED_PENALTY_FUNCTION_LOGCOSHH_CPP: the C++ equivalent to compute_decomposed_penalty_function_logcoshH. --
 * ---------------------------------------------------------------------------------------------------------
 */
void compute_decomposed_penalty_function_logcoshH_cpp(
	Matrix3D *weight_spat,
	Matrix3D *weight_spec,
	Matrix3D *new_estimate,
	Matrix3D *image_last_iter,
	PenaltyParams penalty_params,
	// outputs
	Matrix3D *penalty_function
	)
{
    
	// TODO: is assumption that weight_spat/spec are the same size a valid one?

    penalty_function->zeros(); // initialize
    
    size_t num_rows, num_cols, num_bands;
	weight_spat->getDims(&num_rows, &num_cols, &num_bands);
    
    // Spatial transforms
    if(penalty_params.spat_penalty_flag) {
        for(int row_inc = -1; row_inc <= 1; ++row_inc) {
            for(int col_inc = -1; col_inc <= 1; ++col_inc) {

                // skip the one iter of this we don't care about
                if ( row_inc == 0 && col_inc == 0 ) {
                    continue;
                }

                double neighbor_weight = 1 / sqrt( (double)row_inc * row_inc + col_inc * col_inc);

                // for each spatial point, we consider (r,c,:) of new_est as the curr. pixel and (r-row_inc,c-col_inc,:) of image_last_iter as the neighbor
                #pragma omp parallel for
                for (int r = 0; r < num_rows; ++r) {
                    int r_neigh = r + row_inc;
                    if (r_neigh < 0 || r_neigh >= num_rows) {
                        continue; // neighbor is in row w/o correspondance to any row in image
                    }

                    for (int c = 0; c < num_cols; ++c) {

                        int c_neigh = c + col_inc;
                        if (c_neigh < 0 || c_neigh >= num_cols) {
                            continue; // neighbor is in col w/o correspondance to any col in image
                        }
                        // now we can stop worrying about this point doing any improper indexing

                        // If neighbor_mask would be false, bail on this spatial point

                        double weight_spat_ij = (*weight_spat)(r, c, 0); // for spatial work, all bands at a point will be valid if any are
                        double weight_spat_neigh_ij = (*weight_spat)(r_neigh, c_neigh, 0); // for spatial work, all bands at a point will be valid if any are

                        // make sure we don't increment penalty_* for any point for which neighbour_mask would be false
                        if (weight_spat_ij == 0 || weight_spat_ij != weight_spat_ij) { // if this cell zero or NaN, don't bother
                            continue;
                        }

                        if (weight_spat_neigh_ij == 0 || weight_spat_neigh_ij != weight_spat_neigh_ij) { // if neighbor zero or NaN, don't bother
                            continue;
                        }

                        for (int b = 0; b < num_bands; ++b) {
                            // Increment the pixel corresonding to this (r,c,b) point's penalty_* values as appropriate

                            // bracketterm = 2*new_estimate_padded - c_last_iter_padded - temp;
                            double bracketterm_ijk = 2 * (*new_estimate)(r, c, b) - (*image_last_iter)(r, c, b) - (*image_last_iter)(r_neigh, c_neigh, b);
                            double trigArg_ijk = bracketterm_ijk / penalty_params.delta_spat;

							// scalar = (weight_spat + temp_sen(2:c_num_rows+1,2:c_num_cols+1,:))/2;
							double scalar_ijk = ((*weight_spat)(r, c, b) + (*weight_spat)(r_neigh, c_neigh, b)) / 2.0;

                            double cosh_ijk = cosh(trigArg_ijk);

                            // updating penalty_function
                            if (trigArg_ijk > 500) { // approximate log-cosh by abs value
                                (*penalty_function)(r, c, b) += scalar_ijk * neighbor_weight * penalty_params.beta_spat * penalty_params.delta_spat * penalty_params.delta_spat * trigArg_ijk;
                            }
                            else { // use actual log-cosh
                                (*penalty_function)(r, c, b) += scalar_ijk * neighbor_weight * penalty_params.beta_spat * penalty_params.delta_spat * penalty_params.delta_spat * log(cosh_ijk);
                            }
                        }
                    }
                }

            }
        }
    }

    // Spectral transforms
    if(penalty_params.spec_penalty_flag) {
        int span = penalty_params.spectral_span;
        for (int band_inc = -span; band_inc <= span; ++band_inc) {
            // skip the null transformation
            if (band_inc == 0) {
                continue;
            }

            double neighbor_weight = 1 / (double)abs(band_inc);

            #pragma omp parallel for
            for (int b = 0; b < num_bands; ++b) {

                int b_neigh = b + band_inc;
                if (b_neigh < 0 || b_neigh >= num_bands) {
                    continue; // neighbor is in row w/o correspondance to any row in image
                }

                for (int r = 0; r < num_rows; ++r) {
                    for (int c = 0; c < num_cols; ++c) {
                        
                        // If neighbor_mask would be false, bail on this spatial point

                        double weight_spec_ij = (*weight_spec)(r, c, 0);

                        // make sure we don't increment penalty_* for any spatial location where the sensitivity image would be zero
                        if (weight_spec_ij == 0 || weight_spec_ij != weight_spec_ij) {
                            continue;
                        }

                        // Increment the pixel corresonding to this (r,c,b) point's penalty_* values as appropriate

                        // bracketterm = 2*new_estimate_padded - c_last_iter_padded - temp;
                        double bracketterm_ijk = 2 * (*new_estimate)(r, c, b) - (*image_last_iter)(r, c, b) - (*image_last_iter)(r, c, b_neigh);
                        double trigArg_ijk = bracketterm_ijk / penalty_params.delta_spec;

						double scalar_ijk = (*weight_spec)(r, c, b) + (*weight_spec)(r, c, b_neigh);

                        double cosh_ijk = cosh(trigArg_ijk);

                        // updating penalty_function
                        if (trigArg_ijk > 500) { // approximate log-cosh by abs value
                            (*penalty_function)(r, c, b) += 0.5 * scalar_ijk * neighbor_weight * penalty_params.beta_spec * penalty_params.delta_spec * penalty_params.delta_spec * trigArg_ijk;
                        }
                        else { // use actual log-cosh
                            (*penalty_function)(r, c, b) += 0.5 * scalar_ijk * neighbor_weight * penalty_params.beta_spec * penalty_params.delta_spec * penalty_params.delta_spec * log(cosh_ijk);
                        }
                    }
                }
            }
        }
    }
}

/*
 * -------------------------------------------------------------------------------------------------
 * -- COMPUTE_LOGCOSH_PENALTY_H_CPP: the C++ equivalent to compute_logcosh_penalty_H. --
 * -------------------------------------------------------------------------------------------------
 */
void compute_logcosh_penalty_H_cpp(
	Matrix3D *sensitivity_image, 
	Matrix1D *spec_pnt, // length: num_bands
	Matrix3D *update_factor, 
	Matrix3D *image_last_iter, 
	PenaltyParams penalty_params,
	// outputs
	double *o_penalty,
	double *d_penalty,
	Matrix3D *new_estimate
	)
{

    // NOTE: parent responsible for allocating the outputs
    // (new_estimate is the same size as the other Matrix3D instances)
    
	// if there are negative values, die
	size_t numEls = image_last_iter->num_els();
	for (int i = 0; i < numEls; ++i) {
		if ( (*image_last_iter)(i) < 0) {
			printf("Negative pixel values found in image_last_iter, please check input\n");
			return; // TODO: maybe having a return value is better than this thing Ke did (which lets it die violently later on)
		}
	}

	size_t num_rows, num_cols, num_bands;
	image_last_iter->getDims(&num_rows, &num_cols, &num_bands);

	// For now, let's try to do without image_mask

	// em_estimate_scaled = update_factor.*image_last_iter
	Matrix3D em_estimate_scaled(num_rows, num_cols, num_bands);
    #pragma omp parallel for
	for (int i = 0; i < numEls; ++i) {
		em_estimate_scaled(i) = (*update_factor)(i) * (*image_last_iter)(i);
	}

	if ( !(penalty_params.spat_penalty_flag || penalty_params.spec_penalty_flag) ) {
		// new_image = em_estimate_scaled./sensitivity_image
        #pragma omp parallel for
		for (int i = 0; i < numEls; ++i) {
			(*new_estimate)(i) = em_estimate_scaled(i) / (*sensitivity_image)(i);
		}

		*d_penalty = *o_penalty = 0;
		return;
	}

	// new_estimate = image_last_iter
    #pragma omp parallel for
	for (int i = 0; i < numEls; ++i) {
		(*new_estimate)(i) = (*image_last_iter)(i);
	}

	Matrix3D* weight_spat = sensitivity_image; // create an alias, for convenience
	Matrix3D weight_spec(num_rows, num_cols, num_bands);
	// weight_spec = repmat(reshape(spec_pnt,[1,1,num_bands]),[c_num_rows,c_num_cols,1]);
	for (int r = 0; r < num_rows; ++r) {
		for (int c = 0; c < num_cols; ++c) {
			for (int b = 0; b < num_bands; ++b) {
				weight_spec.ind(r, c, b) = spec_pnt->ind(b) * sensitivity_image->ind(r,c,b);
			}
		}
	}

	double h_val = 0.1;
	Matrix3D h_val_matrix(num_rows, num_cols, num_bands);
	h_val_matrix.fill(h_val);
	
	// saved matrices we'll use through the loop
	Matrix3D penalty_function(num_rows, num_cols, num_bands);
	Matrix3D penalty_gradient(num_rows, num_cols, num_bands);
	Matrix3D penalty_hessian(num_rows, num_cols, num_bands);
	Matrix3D obj_function(num_rows, num_cols, num_bands);
	Matrix3D obj_gradient(num_rows, num_cols, num_bands);
	Matrix3D obj_hessian(num_rows, num_cols, num_bands);
	Matrix3D delta(num_rows, num_cols, num_bands);
	Matrix3D rsm_new_estimate(num_rows, num_cols, num_bands);
	Matrix3D new_penalty_function(num_rows, num_cols, num_bands);

	Matrix3D_bool hit_boundary(num_rows, num_cols, num_bands);
	Matrix3D_bool within_restricted_region(num_rows, num_cols, num_bands);

	/* *********************
	 * *** THE MAIN LOOP ***
	 * *********************/
	for(int iter=0;iter<penalty_params.max_iter;++iter) {
        printf("Restricted Step Method, Iter %d\n",iter+1);
		
		// force a nonnegativity constraint on the restricted region so it ensures that the image
		// update gets nonnegative value

		// TODO: is it just me, or is this crazy? Why does setting h_val_matrix to a lower value enforce non-neg?
		
		// temp = new_estimate - h_val_matrix;
        // h_val_matrix(temp<0) = 0.99*new_estimate(temp<0);
        #pragma omp parallel for
		for (int i = 0; i < numEls; ++i) {
			if( (*new_estimate)(i) - h_val_matrix(i) < 0 ) {
				h_val_matrix(i) = 0.99 * (*new_estimate)(i);
			}
		}
		
		// initializing penalty_*
		penalty_function.zeros();
		penalty_gradient.zeros();
		penalty_hessian.zeros();
		
		// Spatial transforms
		if(penalty_params.spat_penalty_flag) {
			for(int row_inc = -1; row_inc <= 1; ++row_inc) {
				for(int col_inc = -1; col_inc <= 1; ++col_inc) {
					
					// skip the one iter of this we don't care about
					if ( row_inc == 0 && col_inc == 0 ) {
						continue;
					}
					
					double neighbor_weight = 1 / sqrt( (double)row_inc * row_inc + col_inc * col_inc);

					// for each spatial point, we consider (r,c,:) of new_est as the curr. pixel and (r-row_inc,c-col_inc,:) of image_last_iter as the neighbor
					#pragma omp parallel for
                    for (int r = 0; r < num_rows; ++r) {
						int r_neigh = r + row_inc;
						if (r_neigh < 0 || r_neigh >= num_rows) {
							continue; // neighbor is in row w/o correspondance to any row in image
						}

						for (int c = 0; c < num_cols; ++c) {

							int c_neigh = c + col_inc;
							if (c_neigh < 0 || c_neigh >= num_cols) {
								continue; // neighbor is in col w/o correspondance to any col in image
							}
							// now we can stop worrying about this point doing any improper indexing

							// If neighbor_mask would be false, bail on this spatial point

							double sensitivity_ij = (*sensitivity_image)(r, c, 0); // for spatial work, all bands at a point will be valid if any are
							double sensitivity_neigh_ij = (*sensitivity_image)(r_neigh, c_neigh, 0); // for spatial work, all bands at a point will be valid if any are

							// make sure we don't increment penalty_* for any point for which neighbour_mask would be false
							if (sensitivity_ij == 0 || sensitivity_ij != sensitivity_ij) { // if this cell zero or NaN, don't bother
								continue;
							}
                            
                            if (sensitivity_neigh_ij == 0 || sensitivity_neigh_ij != sensitivity_neigh_ij) { // if neighbor zero or NaN, don't bother
                                continue;
                            }

							for (int b = 0; b < num_bands; ++b) {
								// Increment the pixel corresonding to this (r,c,b) point's penalty_* values as appropriate

								// bracketterm = 2*new_estimate_padded - c_last_iter_padded - temp;
								double bracketterm_ijk = 2 * (*new_estimate)(r, c, b) - (*image_last_iter)(r, c, b) - (*image_last_iter)(r_neigh, c_neigh, b);
								double trigArg_ijk = bracketterm_ijk / penalty_params.delta_spat;

								// scalar = (weight_spat + temp_sen(2:c_num_rows+1,2:c_num_cols+1,:))/2;
								double scalar_ijk = ((*weight_spat)(r, c, b) + (*weight_spat)(r_neigh, c_neigh, b)) / 2.0;

								double cosh_ijk = cosh(trigArg_ijk);

								// updating penalty_function
								if (trigArg_ijk > 500) { // approximate log-cosh by abs value
									penalty_function(r, c, b) += scalar_ijk * neighbor_weight * penalty_params.beta_spat * penalty_params.delta_spat * penalty_params.delta_spat * trigArg_ijk;
								}
								else { // use actual log-cosh
									penalty_function(r, c, b) += scalar_ijk * neighbor_weight * penalty_params.beta_spat * penalty_params.delta_spat * penalty_params.delta_spat * log(cosh_ijk);
								}

								// updating penalty_gradient
								penalty_gradient(r, c, b) += 2 * scalar_ijk * neighbor_weight * penalty_params.beta_spat * penalty_params.delta_spat * tanh(trigArg_ijk);

								// updating penalty_hessian
								penalty_hessian(r, c, b) += 4 * scalar_ijk * neighbor_weight * penalty_params.beta_spat / (cosh_ijk * cosh_ijk);
							}

						}
					}
					
				}
			}
		}
		
		// Spectral transforms
		if(penalty_params.spec_penalty_flag) {
			int span = penalty_params.spectral_span;
			for (int band_inc = -span; band_inc <= span; ++band_inc) {
				// skip the null transformation
				if (band_inc == 0) {
					continue;
				}

				double neighbor_weight = 1 / (double)abs(band_inc);
                
                #pragma omp parallel for
				for (int b = 0; b < num_bands; ++b) {

					int b_neigh = b + band_inc;
					if (b_neigh < 0 || b_neigh >= num_bands) {
						continue; // neighbor is in row w/o correspondance to any row in image
					}

					for (int r = 0; r < num_rows; ++r) {
						for (int c = 0; c < num_cols; ++c) {
                            
                            // If neighbor_mask would be false, bail on this spatial point

                            double sensitivity_ij = (*sensitivity_image)(r, c, 0);

                            // make sure we don't increment penalty_* for any spatial location where the sensitivity image would be zero
                            if (sensitivity_ij == 0 || sensitivity_ij != sensitivity_ij) {
                                continue;
                            }

							// Increment the pixel corresonding to this (r,c,b) point's penalty_* values as appropriate

							// bracketterm = 2*new_estimate_padded - c_last_iter_padded - temp;
							double bracketterm_ijk = 2 * (*new_estimate)(r, c, b) - (*image_last_iter)(r, c, b) - (*image_last_iter)(r, c, b_neigh);
							double trigArg_ijk = bracketterm_ijk / penalty_params.delta_spec;

							double scalar_ijk = weight_spec.ind(r, c, b) + weight_spec.ind(r, c, b_neigh);

							double cosh_ijk = cosh(trigArg_ijk);

							// updating penalty_function
							if (trigArg_ijk > 500) { // approximate log-cosh by abs value
								penalty_function(r, c, b) += 0.5 * scalar_ijk * neighbor_weight * penalty_params.beta_spec * penalty_params.delta_spec * penalty_params.delta_spec * trigArg_ijk;
							}
							else { // use actual log-cosh
								penalty_function(r, c, b) += 0.5 * scalar_ijk * neighbor_weight * penalty_params.beta_spec * penalty_params.delta_spec * penalty_params.delta_spec * log(cosh_ijk);
							}

							// updating penalty_gradient
							penalty_gradient(r, c, b) += scalar_ijk * neighbor_weight * penalty_params.beta_spec * penalty_params.delta_spec * tanh(trigArg_ijk);

							// updating penalty_hessian
							penalty_hessian(r, c, b) += 2 * scalar_ijk * neighbor_weight * penalty_params.beta_spec / (cosh_ijk * cosh_ijk);


						}
					}
				}
			}
		}
		
		// penalty_function(~image_mask) = 0;
        // penalty_gradient(~image_mask) = 0;
        // penalty_hessian(~image_mask) = 0;
        // obj_function = -em_estimate_scaled.*log(new_estimate) + sensitivity_image.*new_estimate + penalty_function;
        // obj_gradient = -em_estimate_scaled./new_estimate + sensitivity_image + penalty_gradient;
        // obj_heissian = em_estimate_scaled./(new_estimate.^2) + penalty_hessian;
		// delta_raw = -obj_gradient./obj_heissian;%heissian matrix is actually a diagonal matrix so inversion becomes element-wise division
        // delta = delta_raw;
		// %evaluate the obj_function at new estimate and compute r
		// rsm_new_estimate = new_estimate + delta;

		// NOTE: we won't actually zero out parts of penalty_* since we never need them again
		
        #pragma omp parallel for
		for (int i = 0; i < numEls; ++i) {
			double penalty_function_zeroed_ijk = ((*sensitivity_image)(i) > 0) ? penalty_function(i) : 0; // use zero if image_mask is false here
			double penalty_gradient_zeroed_ijk = ((*sensitivity_image)(i) > 0) ? penalty_gradient(i) : 0;
			double penalty_heissian_zeroed_ijk = ((*sensitivity_image)(i) > 0) ? penalty_hessian(i) : 0;
			
			obj_function(i) = (-em_estimate_scaled(i) * log((*new_estimate)(i))) + \
				((*sensitivity_image)(i) * (*new_estimate)(i)) + penalty_function_zeroed_ijk;
			obj_gradient(i) = (-em_estimate_scaled(i) / (*new_estimate)(i)) + \
				(*sensitivity_image)(i) + penalty_gradient_zeroed_ijk;
			obj_hessian(i) = em_estimate_scaled(i) / ( (*new_estimate)(i) * (*new_estimate)(i) ) + \
				penalty_heissian_zeroed_ijk;
			
			delta(i) = -obj_gradient(i) / obj_hessian(i);

			// use the raw version of delta to populate our boolean matrices, hit_boundary and within_restricted_region
			within_restricted_region(i) = abs(delta(i)) <= h_val_matrix(i);
			hit_boundary(i) = !within_restricted_region(i);

			// updating this pixel of delta
			if (delta(i) > h_val_matrix(i)) {
				delta(i) = h_val_matrix(i);
			}
			else if (delta(i) < -h_val_matrix(i)) {
				delta(i) = -h_val_matrix(i);
			}

			rsm_new_estimate(i) = (*new_estimate)(i) + delta(i);
		}

		// get decomposed penalty value
		// we are calling the return value new_penalty_function, and when we apply the image mask, we will still call it this (instead of temp)

		// temp = image_mask.*compute_decomposed_penalty_function_logcoshH(weight_spat,weight_spec, rsm_new_estimate, image_last_iter,beta_spat, delta_spat, beta_spec, delta_spec,spectral_span, spat_penalty_on, spec_penalty_on);
		//obj_function_quadratic_approximate = obj_function + obj_gradient.*delta + 0.5*obj_heissian.*delta.*delta;
		//obj_function_2 = -em_estimate_scaled.*log(rsm_new_estimate) + sensitivity_image.*rsm_new_estimate + temp;
		//%compute r
		//obj_diff = obj_function - obj_function_2;
		//obj_approximate_diff = obj_function - obj_function_quadratic_approximate;
		//r_vals = obj_diff. / obj_approximate_diff;
		//%if r_vals>0.75, then the quadratic approximates the function well.
		//%and if the delta is within the restricted region, then the h_val is good enough
		//h_val_matrix(r_vals>0.75 & within_restricted_region) = h_val_matrix(r_vals>0.75 & within_restricted_region);
		//%if the delta is at the boundary of the region, then we need to
		//%increase the region
		//h_val_matrix(r_vals>0.75 & hit_boundary) = 2 * h_val_matrix(r_vals>0.75 & hit_boundary);
		//%if r_vals<0.25, bad approximation, the curvature is too flat, we
		//%need to use a smaller restricted region
		//h_val_matrix(r_vals<0.25) = 0.25*abs(delta(r_vals<0.25));
		//new_estimate(r_vals>0) = rsm_new_estimate(r_vals>0);
        
		compute_decomposed_penalty_function_logcoshH_cpp(weight_spat, &weight_spec, &rsm_new_estimate, image_last_iter, penalty_params, &new_penalty_function);

		double objective_function_quadratic_approximate_ijk, obj_function_2_ijk, obj_diff_ijk, obj_approximate_diff_ijk, r_vals_ijk;
		
        #pragma omp parallel for
        for (int i = 0; i < numEls; ++i) {
			new_penalty_function(i) = ((*sensitivity_image)(i) > 0) ? new_penalty_function(i) : 0;

			objective_function_quadratic_approximate_ijk = obj_function(i) + obj_gradient(i)*delta(i) + 0.5*obj_hessian(i)*delta(i)*delta(i);
			obj_function_2_ijk = -em_estimate_scaled(i) * log(rsm_new_estimate(i)) + (*sensitivity_image)(i) * rsm_new_estimate(i) + new_penalty_function(i);
			
			obj_diff_ijk = obj_function(i) - obj_function_2_ijk;
			obj_approximate_diff_ijk = obj_function(i) - objective_function_quadratic_approximate_ijk;

			r_vals_ijk = obj_diff_ijk / obj_approximate_diff_ijk;

			// modifying our variables that persist over iterations: h_val_matrix and new_estimate

			// there was an r_val > 0.75 and within_restricted_region case in Matlab, but it did a no-op
			if (r_vals_ijk > 0.75 && hit_boundary(i)) {
				h_val_matrix(i) = 2 * h_val_matrix(i);
			}
			else if (r_vals_ijk < 0.25) {
				h_val_matrix(i) = 0.25 * abs(delta(i));
			}

			if (r_vals_ijk > 0) {
				(*new_estimate)(i) = rsm_new_estimate(i);
			}
		}
        
        printf("Trust regions & Image Estimate updated\n");

	}
    
    // the new estimate is in the output as expected

    // produce d_penalty and o_penalty
	{ // a scope to remove temporary matrix as soon as possible
		Matrix3D new_penalty_value(num_rows, num_cols, num_bands);

		compute_decomposed_penalty_function_logcoshH_cpp(weight_spat, &weight_spec, new_estimate, image_last_iter, penalty_params, &new_penalty_value);
		*d_penalty = 0;
		for (int i = 0; i < numEls; ++i) {
			(*d_penalty) += new_penalty_value(i);
		}
	}

	{ // a scope to remove temporary matrices as soon as possible
		Matrix3D new_penalty_value_spat(num_rows, num_cols, num_bands);
		Matrix3D new_penalty_value_spec(num_rows, num_cols, num_bands);
		compute_original_penalty_function_spat_cpp(sensitivity_image, new_estimate, penalty_params, &new_penalty_value_spat);
		compute_original_penalty_function_spec_cpp(sensitivity_image, new_estimate, penalty_params, &new_penalty_value_spec);
		*o_penalty = 0;
		for (int i = 0; i < numEls; ++i) {
			(*o_penalty) += (*weight_spat)(i)*new_penalty_value_spat(i);
			(*o_penalty) += weight_spec(i)*new_penalty_value_spec(i);
		}
	}

}

/* 
 * A MEX version of compute_new_image_logcosh_prior.m, build for both
 * time & memory efficiency.
 *
 * Inputs:
 * 0. sensitivity_image
 * 1. update_factor
 * 2. image_last_iter
 * 3. beta_spat
 * 4. delta_spat
 * 5. beta_spec
 * 6. delta_spec
 * 7. spectral_span
 * 8. max_iter
 * 9. spat_penalty_on
 * 10. spec_penalty_on
 *
 * Outputs:
 * 0. o_penalty
 * 1. d_penalty
 * 2. new_image
 *
 */

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    #define sensitivity_image_in    prhs[0]
	#define spec_pnt_in				prhs[1]
    #define update_factor_in        prhs[2]
    #define image_last_iter_in      prhs[3]
    #define beta_spat_in            prhs[4]
    #define delta_spat_in           prhs[5]
    #define beta_spec_in            prhs[6]
    #define delta_spec_in           prhs[7]
    #define spectral_span_in        prhs[8]
    #define max_iter_in             prhs[9]
    #define spat_penalty_on_in      prhs[10]
    #define spec_penalty_on_in      prhs[11]
	#define omp_threads_in          prhs[12]
    
    #define o_penalty_out           plhs[0]
    #define d_penalty_out           plhs[1]
    #define new_image_out           plhs[2]
    
    // Checking the properties of the inputs & outputs
    // TODO?

    // Gathering input values into form the c++ function expects
    
    
    mwSize numDims = mxGetNumberOfDimensions(image_last_iter_in);
    // Note: the docs for R2016a suggest that the size_t here should be mwSize, but in Matlab R2015a, it is size_t. Perhaps not a solid API.
    const size_t* dims = mxGetDimensions(image_last_iter_in);
    if(numDims != 3) {
        printf("We've got a problem...\n");
        return;
    }
    size_t c_num_rows = dims[0];
    size_t c_num_cols = dims[1];
    size_t num_bands = dims[2];
    
    PenaltyParams penalty_params;
    
    penalty_params.beta_spat = mxGetScalar(beta_spat_in);
    penalty_params.delta_spat = mxGetScalar(delta_spat_in);
    penalty_params.beta_spec = mxGetScalar(beta_spec_in);
    penalty_params.delta_spec = mxGetScalar(delta_spec_in);
    penalty_params.spat_penalty_flag = mxGetScalar(spat_penalty_on_in)==1; // conversion from double -> bool
    penalty_params.spec_penalty_flag = mxGetScalar(spec_penalty_on_in)==1; // conversion from double -> bool
    penalty_params.max_iter = (int)mxGetScalar(max_iter_in);
    penalty_params.spectral_span = (int)mxGetScalar(spectral_span_in);
    
    // Call the c++ function
    
    double o_penalty, d_penalty;
    
    // Wrap the input matrices in the Matrix3D class (Without allocating space; we already have them in memory)
    Matrix3D sensitivity_image(c_num_rows,c_num_cols,num_bands,false);
    sensitivity_image.setData((double*)mxGetData(sensitivity_image_in));
    sensitivity_image.setDataForfeit(false); // Don't destroy the input data

	Matrix1D spec_pnt(num_bands, false);
	spec_pnt.setData((double*)mxGetData(spec_pnt_in));
	spec_pnt.setDataForfeit(false); // Don't destroy the input data
    
    Matrix3D update_factor(c_num_rows,c_num_cols,num_bands,false);
    update_factor.setData((double*)mxGetData(update_factor_in));
    update_factor.setDataForfeit(false); // Don't destroy the input data
    
    Matrix3D image_last_iter(c_num_rows,c_num_cols,num_bands,false);
    image_last_iter.setData((double*)mxGetData(image_last_iter_in));
    image_last_iter.setDataForfeit(false); // Don't destroy the input data
    
    // Build a shell around a Matlab-allocated location for the output matrix
	new_image_out = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // same size & shape as the inputs
    Matrix3D new_estimate(c_num_rows,c_num_cols,num_bands,false);
    new_estimate.setData(mxGetPr(new_image_out)); // use Matlab's manager to alloc so we can return it from this function
    new_estimate.setDataForfeit(false); // allows us to carry out our data
    
	// setup OpenMP
	int omp_threads = (int)std::floor(mxGetScalar(omp_threads_in));
	setupOpenMP(omp_threads);

	// Call the c++ function
    compute_logcosh_penalty_H_cpp(&sensitivity_image, &spec_pnt, &update_factor, &image_last_iter, penalty_params, &o_penalty, &d_penalty, &new_estimate);
    
    // Set outputs for MATLAB
    o_penalty_out = mxCreateDoubleScalar(o_penalty);
    d_penalty_out = mxCreateDoubleScalar(d_penalty);
    
}
