function f = calculate_spatial_kernel_anygeo(geoflags,in_matrix,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,pix_spacing, mro_position,pix_positions,pix_normals, use_mex)
% Update by Linyun 03/02/2016
% Using huge kernel as variable
tic;

% Run MEX version if we have it
if use_mex && ~geoflags.precompute
    if exist('general_proj_spat_mex','file')
        global mex_num_threads;
        if(isempty(mex_num_threads))
            mex_num_threads=-1; % default to max threads
        end
        
        f = general_proj_spat_mex(in_matrix,geoflags.forward,geoflags.offnadir,ddr_iof,ddr_lat,ddr_lon,sb,size_kernel,lat_number,lon_number,pix_spacing,mro_position,pix_positions,pix_normals,mex_num_threads);
        toc;
        return;
    else
        disp('C++ version of spatial projections not found for platform. Using Matlab versions...');
    end
end

padding = (size_kernel-1)/2;
gaussian_kernel = ones(size_kernel,size_kernel,num_bands);
num_elements = size_kernel * size_kernel;

if geoflags.precompute
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %Initialize the kernel to write to file. The order of indexing is
    %row_iter,col_iter,band,row,col
    huge_kernel = zeros(size_kernel, size_kernel, num_bands, num_rows, num_cols); 
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
else
    
    if geoflags.forward
        fragment = ones(size_kernel,size_kernel,num_bands);
        conv_layer = ones(num_rows, num_cols, num_bands);
    else
        factor = ones(c_num_rows, c_num_cols, num_bands); %cube to hold factor for backprojecting
        conv_layer = zeros(c_num_rows+2*padding, c_num_cols+2*padding, num_bands);
    end
end
%-------------------physical constants-----------------------------
ifov = 61.5e-6;
r = 3396190.0;
spectral_sampling = 6.55; %nm/channel
nominal_altitude_sq = (266e3)^2; %nominal altitude of MRO above the surface
%-------------------projected areas-----------------------------
projected_areas = zeros(1,c_num_rows);
if geoflags.offnadir
    %calculate projected area (assuming spherical planet and radius one) of
    %each row of pixels, given that the grid is latitude- and longitude-regular.
    %Will be adjusted on a per-pixel basis to account
    %for local variations of elevation, slope and viewing angle.
    lon_spacing = (ddr_lon(c_num_cols)-ddr_lon(1))/(c_num_cols-1);
    lat_spacing = (ddr_lat(c_num_rows)-ddr_lat(1))/(c_num_rows-1);
    for c_row = 1:c_num_rows
        top_lat = ddr_lat(c_row) - 1/2*lat_spacing;
        bottom_lat = ddr_lat(c_row) + 1/2*lat_spacing;
        projected_areas(c_row) = abs(lon_spacing*(cos(pi/2-top_lat) - cos(pi/2-bottom_lat)));
    end
else
    projected_areas(:,:) = pix_spacing^2; %in old code this was multiplied by 4/(4*pi), but the 4s cancel and the 1/pi is exported to the normalization coefficient (leading_coeff_default)
end

for col = 1:num_cols
    fwhm_alt1_sq = (sb(1,col,1:num_bands)/spectral_sampling*ifov).^2; %the fwhm squared if altitude = 1 (will rescale)        
    for row = 1:num_rows
        c_row = lat_number(row,col); 
        c_col = lon_number(row,col);
        apparent_areas = zeros(size_kernel,size_kernel); %Will multiply the kernel by this array, elementwise, before normalizing
        %scale fwhm according to distance from MRO to target
        if geoflags.offnadir
            R_target = reshape(pix_positions(c_row,c_col,:), [1,3]); %position vector of the target pixel
            LOS_target = R_target - mro_position(row,:); %Line of sight from MRO to the target pixel
            %Compute some default values for out of bounds pixels 
            normal_target = reshape(pix_normals(c_row,c_col,:), [1,3]); %normal at the surface from the target
            viewing_area_cosine_target = abs(LOS_target*normal_target'/sqrt(LOS_target*LOS_target')); %normal_target is already normalized. This is the same as cosine of emission angle
            local_radius_sq = R_target*R_target';
            slope_cosine_target = abs(R_target*normal_target'/sqrt(R_target*R_target')); %normal_target is already normalized
            apparent_areas(padding+1,padding+1) = local_radius_sq*projected_areas(c_row)*viewing_area_cosine_target/slope_cosine_target;
            LOS_target_sq = LOS_target*LOS_target';
        else
            LOS_target_sq = nominal_altitude_sq;
            pix_lat = ddr_iof(row,col,4)*pi/180;
            pix_lon = ddr_iof(row,col,5)*pi/180;
        end
        fwhm_sq = fwhm_alt1_sq*LOS_target_sq;
        leading_coeff = log(16)./(pi*fwhm_sq);
        leading_coeff_big = repmat(leading_coeff,[size_kernel,size_kernel,1]);
        
        if ~geoflags.precompute && geoflags.forward
            fragment = safeSubset(in_matrix, c_row-padding,c_row+padding, c_col-padding,c_col+padding, 0);
        end
        
        % Fully vectorized kernel construction, nadir and off-nadir
        if ~geoflags.offnadir
            iter_lats = safeSubset1d(ddr_lat,c_row-padding,c_row+padding,nan);
            iter_lats = repmat(iter_lats,[size_kernel,1])';
            
            iter_lons = safeSubset1d(ddr_lon,c_col-padding,c_col+padding,nan);
            iter_lons = repmat(iter_lons,[size_kernel,1]);
            
            ys = r*(pix_lat-iter_lats);
            r_effs = r*sin(pi/2-abs(pix_lat+iter_lats)/2);
            xs = r_effs.*(pix_lon-iter_lons);
            
            
            
            dist_sqs = xs.^2+ys.^2;
            
            dist_sqs_big = repmat(dist_sqs,[1,1,num_bands]);
            
            gaussian_kernel = leading_coeff_big.*exp(-(pi*leading_coeff_big) .* dist_sqs_big);
            
            out_of_bounds_positions = isnan(gaussian_kernel);
            gaussian_kernel(out_of_bounds_positions) = pi * leading_coeff_big(out_of_bounds_positions);
            
            apparent_areas(:,:) = pix_spacing^2; % == projected_areas in off-nadir case
            apparent_areas(isnan(xs)|isnan(ys)) = 1; %This makesthe out-of-bounds values much lower
            
        else % geoflags.offnadir is true
            R_neighbors = safeSubset(pix_positions,c_row-padding,c_row+padding,c_col-padding,c_col+padding,NaN);
            LOS_neighbors = R_neighbors - repmat(permute(mro_position(row,:),[1 3 2]),[size_kernel, size_kernel, 1]);
            LOS_neighbor_sqs = dot_fast(LOS_neighbors,LOS_neighbors,3);
            %find distance between target and projection of
            %neighbor onto plane intersecting target and
            %perpendicular to LOS_target
            LOS_target_big = repmat(permute(LOS_target,[1 3 2]),[size_kernel, size_kernel, 1]);
            
            cossqs = dot_fast(LOS_target_big,LOS_neighbors,3).^2 ./ (LOS_target_sq*LOS_neighbor_sqs); %cosine squared of the angle between the LOS's
            tansqs = 1./cossqs-1;
            ttn_dist_sqs = LOS_target_sq*tansqs;
            %compute the (normalized) gaussian:
            gaussian_kernel = leading_coeff_big .*exp(-(pi*leading_coeff_big) .* repmat(ttn_dist_sqs,[1,1,num_bands]));
            
            % filling in out-of-bounds areas
            out_of_bounds_positions = isnan(gaussian_kernel);
            gaussian_kernel(out_of_bounds_positions) = pi * leading_coeff_big(out_of_bounds_positions);
            
            %populate apparent area matrix
            normal_neighbors = safeSubset(pix_normals,c_row-padding,c_row+padding,c_col-padding,c_col+padding,NaN);
            
            viewing_area_cosines = abs(dot_fast(LOS_neighbors,normal_neighbors,3)./sqrt(LOS_neighbor_sqs)); %normal_neighbor is already normalized
            local_radius_sqs = dot_fast(R_neighbors,R_neighbors,3);
            slope_cosines = abs(dot_fast(R_neighbors,normal_neighbors,3)./sqrt(local_radius_sqs)); %normal_neighbor is already normalized
            projected_areas_subset = safeSubset1d(projected_areas,c_row-padding,c_row+padding,NaN)';
            apparent_areas = local_radius_sqs.*repmat(projected_areas_subset,[1,size_kernel]).*viewing_area_cosines./slope_cosines;
            
            % filling in out-of-bounds areas
            out_of_bounds_positions = isnan(apparent_areas);
            apparent_areas(out_of_bounds_positions) = 1;
        end
        %multiply Gaussian kernel by the apparent area matrix
        gaussian_kernel = gaussian_kernel.*repmat(apparent_areas, [1,1,num_bands]);
        if geoflags.offnadir
            %Normalize Gaussian kernel
            gaussian_kernel_vol = sum(sum(gaussian_kernel(:,:,:))); %a 1x1x246 array, each entry the volume of the corresponding band
            gaussian_kernel = gaussian_kernel./repmat(gaussian_kernel_vol,[size_kernel,size_kernel,1]);
        end
        if geoflags.precompute
            %*****************Compute total kernel**********************************
            %Convolve the two kernels (gaussian_kernel and smear_kernel) together to get the total PSF
            huge_kernel(:, :, :, row, col) = gaussian_kernel;
            %*****************Total kernel computed!***********************
        else
            if geoflags.forward
                kernel_temp = reshape(gaussian_kernel, num_elements, num_bands);
                fragment_temp = reshape(fragment, num_elements, num_bands);
                conv_layer(row,col,:) = dot_fast(kernel_temp,fragment_temp);
            else
                temp2 = repmat (in_matrix(row, col, 1:num_bands), [size_kernel size_kernel 1]);
                
                temp = temp2.*gaussian_kernel;
                
                %boundaries of patch in the new image that the kernel
                %corresponds to
                row_min = lat_number(row,col);
                row_max = lat_number(row,col)+2*padding;
                col_min = lon_number(row,col);
                col_max = lon_number(row,col)+2*padding;
                
                %compute section of scene
                conv_layer(row_min:row_max,col_min:col_max, :) = conv_layer(row_min:row_max,col_min:col_max, :) + temp;
            end
        end
    end
end

if geoflags.precompute
    f = huge_kernel;
else
    if geoflags.forward
        f = conv_layer;
        display('Spatial Forward Projection: ');
    else
        factor(:, :, :) = conv_layer(1+padding:c_num_rows+padding, 1+padding:c_num_cols+padding, :);
        'back_spat done.';
        f = factor;
        display('Spatial Backprojection: ');
    end
end

toc;
end