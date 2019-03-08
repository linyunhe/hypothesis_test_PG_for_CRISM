% Now encapsulated in a function to
%   1) Make it runnable on CHPC cluster,
%   2) Prevent polution of the namespace
%
% It can have no arguments, so it can still be executed by pressing "run."
%
% ARGUMENTS: Takes pairs of arguments of the form
% 'ParameterName','ParameterValue'. The name should be that of some field
% of params found at the top of the function, and the value should be of
% whatever matlab type is appropriate. If arguments are given, they are
% used instead of whatever values are written into the script; the script
% values can be considered the default values.
function [Cost_function] =test5v5_compact(varargin)

format;

tic_total = tic;

%%%%%%%%%%%%%
%%edit this%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: to do an "X iterations off, Y iterations on" style run, leave the 
% penalty flags off  and set the switchover settings so to turn them on later.
params.spatial_penalty_flag = 0; % 1 for using penalty, 0 for no penalty
params.spectral_penalty_flag = 0; % 1 for using penalty, 0 for no penalty
params.beta_spat = 0.01;
params.delta_spat = 4.0;
params.beta_spec = 0.04;
params.delta_spec = 0.9;
params.max_iter = 3; %maximum iterations be used in Restricted Region Newton's method
params.spectral_span = 3; %number of spectral neighbours considered for spectral penalty

params.is_L_data = 1; % true for L, false for S. For choosing spectral transfer method

params.size_kernel = 11; %must be odd
params.pix_spacing = 12.0; %size of pixel diameter during processing and on local-cartesian output
params.pix_spacing_output = 12.0; %size of pixel diameter in map-projected output
params.pix_value = 0.6; %default intial c pixel value. Must be greater than 0 or less than 1.
params.num_iter = 25; %number of times to iterate through the algorithm.

params.g1_flag = true; % if true, use 1 spectral gaussian even for L data, otherwise use 3 assym. gaus. in that case

%------------Spatial TF parameters-------------------------
params.geoflags.offnadir = 0; % If true: use flat, off-nadir geo. If false: use nadir geo.
%The following are only applicable if using new geometry
params.geoflags.flat = 1; %Approximating a flat surface (still accounting for viewing angle)
%-----------------------------------------------------------

% Switchover feature. If true, the given fragment of code is executed when the given iteration finishes
params.switchover_flag = 1;
params.switchover_at_iter = 5; % num of the iteration after which we apply the switchover
params.switchover_code = 'params.spatial_penalty_flag=1;params.spectral_penalty_flag=1;';


if params.spatial_penalty_flag || params.spectral_penalty_flag
    start_iter = 1;
elseif params.switchover_flag
    start_iter = params.switchover_at_iter + 1;
else 
    start_iter = 1;
end

params.spectral_offset = 0; % in units nm

params.save_c_on_iters_flag = 1; % whether to save the scene on the iters listed below
params.save_mu_on_iters_flag = 0; % whether to save sensor space on the iters listed below
params.save_proj_on_iters_flag = 0;
params.save_input_flag = 0;
params.save_spectrogram_flag = 0;

params.saving_iters = [25]; % note: both c and mu are always saved after the last iteration regardless of these settings

% input files
params.crism_iof_filename = ''; %measured radiance
params.ddr_iof_filename = ''; %gives latitude, logitude, emission angle, etc
params.sb_filename = ''; %parameters for spectral TF
params.wa_filename = ''; %bandpass centers

% input files only necessary for off-nadir geometry
params.filenames.radius_filename = '';
params.filenames.topo_filename = '';
params.filenames.spice_filename = '';

% For starting a run from a scene produced in a previous run.
% If you use this, it is best to still give the original IOF/SSA file above, so the
% algorithm converges to an image estimating it, rather than the intermediate product.
params.start_from_scene_flag = 0;
params.start_from_scene_filename = ''; % an ENVI-format scene that we use as initial guess

% Paths to prepend to input and output file paths respectively ("." is the
% current directory). If empty string, assumes an acceptable full path has
% been given.
params.input_path_start_in = '';
params.output_path_start_in = '';

% input file subsetting
% (all indexed from 1)
% if the *_use_all flag is false, uses the *_subset_bounds as the min and max extent
% only datasets with 640 columns are supported (and matching # of rows & bands in SSA & DDR)
params.rows_use_all = 0;
params.rows_subset_bounds = [1,170];
params.cols_use_all = true;
params.cols_subset_bounds = [250,350];
params.bands_bounds = [193,430]; % the bands that corresond to the range found in the SSA
params.pre_end = 36;
params.post_start = 50;
params.spec_range = 10;
% output files for this run will get this as a prefix to their name
params.observation_base = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user-editable parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUARANTINED PARAMETERS: not of interest to end-user & slated for removal
params.geoflags.precompute = false; %Option to pre-compute the spatial PSF kernel and reference it for future iterations, rather than re-compute it

% Parameters that are not user-controlled (not in GUI, etc)
params.max_iter = 3; % max iters in main penalty loop (in compute_ratio_penalty_H_linear_spec)
% TODO: next one might be worth putting in GUI
params.use_mex_versions = true; % Use MEX C++ translations when possible
params.mex_num_threads = -1; % Number of threads to request. -1 requests the maximum number.

global mex_num_threads; mex_num_threads = params.mex_num_threads;

% ddr band names: for use in outputting cut DDRs
ddr_band_names = {'INA at areoid; deg','EMA at areoid; deg','Phase angle; deg',...
    'Latitude; areocentric; deg N','Longitude; areocentric; deg E',...
    'INA at surface from MOLA; deg','EMA at surface from MOLA; deg',...
    'Slope magnitude from MOLA; deg','MOLA slope azimuth; deg clkwise from N',...
    'Elevation; meters relative to MOLA','Thermal inertia; J m^-2 K^-1 s^-0.5',...
    'Bolometic albedo','Local solar time; hours','Spare'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use function args to replace some params %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = setParamsFromArgs( varargin, params );

[version,algo_type] = get_version_number(); % algorithm & version num

fprintf('Beginning file reading.\n');
fprintf('The time is %s',info_string(params.observation_base,version,algo_type));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data loading from the files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% read data files
crism_iof = read_envi_data(fullfile(params.input_path_start_in,params.crism_iof_filename)); %measured radiance
ddr_iof = read_envi_data(fullfile(params.input_path_start_in,params.ddr_iof_filename)); %gives latitude, logitude, emission angle, etc
sb = read_envi_data(params.sb_filename); %parameters for spectral TF
wa = read_envi_data(params.wa_filename); %bandpass centers

if ~params.is_L_data
   try 
       crism_iof = crism_iof(:,:,params.bands_bounds(1):params.bands_bounds(2));
   catch 
        error('Please examine how many bands the S SB & WA files have. Does not match band bounds given.');
   end
end

if params.start_from_scene_flag
    % An ENVI-format scene that we use as initial guess.
    % Currently must be the same size as the scene crism_setup would create
    % for us
    start_from_scene = read_envi_data(fullfile(params.input_path_start_in,params.start_from_scene_filename));
	start_from_scene = transformScene(start_from_scene,'inverse'); % rotate 180 deg to get to working orientation
end

% checking that a few dimesions of the data with fixed lengths are correct
if params.is_L_data
    sb_num_rows = 10;
else
    sb_num_rows = 1;
end

total_cols = 640; % determined by instrument

if size(sb,1)~= sb_num_rows
    error('SB must have 10 rows for L data and 1 row for S data')
elseif size(ddr_iof,3)~=14
    error('DDR must have 14 bands');
elseif size(wa,1) ~= 1
    error('WA must have 1 row');
elseif size(crism_iof,2) ~= total_cols
    error('IOF/SSA must have 640 columns (currently: %i)',size(crism_iof,2));
elseif size(ddr_iof,2) ~= total_cols
    error('DDR must have 640 columns (currently: %i)',size(ddr_iof,2));
elseif size(sb,2) ~= total_cols
    error('SB must have 640 columns (currently: %i)',size(sb,2));
elseif size(wa,2) ~= total_cols
    error('WA must have 640 columns (currently: %i)',size(wa,2));
end

% checking that the sizes of the subsetted data are consistent
if size(crism_iof,1)~=size(ddr_iof,1)
    error('IOF/SSA and DDR must have same number of rows (currently: %i and %i)',size(crism_iof,1),size(ddr_iof,1));
elseif size(sb,3)~=size(wa,3)
    error('SB and WA must have same number of bands (currently: %i and %i)',size(sb,3),size(wa,3));
end

if params.is_L_data
    col_furthest_bounds = [32,631];
else % S data
    col_furthest_bounds = [26,626];
end

% auto-subsetting from parameter values
if params.rows_use_all
    row_min = 1;
    row_max = size(crism_iof,1);
else
    row_min = params.rows_subset_bounds(1);
    row_max = params.rows_subset_bounds(2);
end

if params.cols_use_all % don't take all 640, but all the valid ones anyway
    col_min = col_furthest_bounds(1);
    col_max = col_furthest_bounds(2);
else
    col_min = params.cols_subset_bounds(1);
    col_max = params.cols_subset_bounds(2);
end

band_min = params.bands_bounds(1);
band_max = params.bands_bounds(2);

% Maybe some of the columns are invalid, if an incomplete dataset was
% padded. Get rid of those, too.

validCols = find(squeeze(crism_iof(1,:,1)) ~= 65535);

if col_min < min(validCols)
    col_min = min(validCols);
end

if col_max > max(validCols)
    col_max = max(validCols);
end

% make sure no columns requested that are out-of-bounds
if (col_min < col_furthest_bounds(1)) || (col_furthest_bounds(2) < col_max)
    error('Invalid columns selected [%i, %i]. Must be within (inclusive) range [%i, %i].', ...
        col_min,col_max, ...
        col_furthest_bounds(1), col_furthest_bounds(2));
end

% one last consistency check: is the correct number of bands selected?
if band_max - band_min + 1 ~= size(crism_iof,3)
    error('The selected band range [%i, %i] does not have the correct number of bands (%i; currently: %i)',band_min,band_max,size(crism_iof,3),band_max-band_min+1);
end

% decide which default bands to write into ENVI files (RGB; 1-indexed like
% in WA & SB files)
if params.is_L_data
    default_bands = [206, 361, 429];
else % S-data
    default_bands = [54, 37, 27];
end

% Determine if any edge rows or cols need to be cut (because they are
% uniformly zero)
has_contents_rows = all(any(crism_iof,2),3);
has_contents_cols = all(any(crism_iof),3);

% only consider the part we've selected so far
has_contents_rows = has_contents_rows(row_min:row_max);
has_contents_cols = has_contents_cols(col_min:col_max);

% determine the number of rows/cols on each side that ought to be cut
num_bad_early_rows = find(has_contents_rows==1,1,'first') - 1;
num_bad_late_rows = length(has_contents_rows) - find(has_contents_rows==1,1,'last');
num_bad_early_cols = find(has_contents_cols==1,1,'first') - 1;
num_bad_late_cols = length(has_contents_cols) - find(has_contents_cols==1,1,'last');

% adjust mins and maxes accordingly

if length(num_bad_early_rows) ~= 0
    row_min = row_min + num_bad_early_rows;
end
if length(num_bad_late_rows) ~=0
    row_max = row_max - num_bad_late_rows;
end
if length(num_bad_early_cols) ~= 0
    col_min = col_min + num_bad_early_cols;
end
if length(num_bad_late_cols) ~= 0
    col_max = col_max - num_bad_late_cols;
end
 

crism_iof = crism_iof(row_min:row_max,col_min:col_max,:); % currently, no support for subsetting bands within SSA
ddr_iof = ddr_iof(row_min:row_max,col_min:col_max,:);
sb = sb(:,col_min:col_max,band_min:band_max);
wa = wa(:,col_min:col_max,band_min:band_max);

% capture dimension sizes in variables that the code uses
[num_rows, num_cols, num_bands] = size(crism_iof);

if any(default_bands < band_min | default_bands > band_max)
    % at least one default band is outside our selection, so don't write
    % them into the file
    default_bands = [];
    disp('Warning: default bands for ENVI headers are out of range for selected data and will not be added to files.');
end

% Adjust default_bands to be relative to the bands that will be in
% output files
default_bands = default_bands - (band_min - 1);

% reshaping & correcting wa
wa = reshape(wa,[num_cols,num_bands]); % for convenience, flatten to 2D
wa = wa + params.spectral_offset;

if params.geoflags.offnadir
    % subsetting the same way as we do the rows of ddr_iof
    mro_position = spice_info(params.filenames.spice_filename,row_min,row_max);
else
    mro_position = 0; % dummy value; if we're not using off-nadir method we don't need this variable
end

% CRISM_IOF negative suppression: doesn't error out, but informs user
if any(crism_iof(:)<0) 
    num_bad_cells = sum(crism_iof(:)<0);
    crism_iof(crism_iof<0) = 0;
    fprintf('WARNING: crism_iof contained %i negative elements (%.2f%% of total), which were suppressed automatically. Check your data set.\n',num_bad_cells,100*num_bad_cells/length(crism_iof(:)));
end

bounds.rowBounds = [row_min, row_max];
bounds.colBounds = [col_min, col_max];
bounds.bandBounds = [band_min, band_max];

% build a description we can use in header files
desc_text = generateParamListDesc(params, version, algo_type, bounds);

fprintf('Performing algorithm on rows %i to %i (%i), cols %i to %i (%i), bands %i to %i (%i) (all 1-indexed).\n',row_min,row_max,num_rows,col_min,col_max,num_cols,band_min,band_max,num_bands);

kernel_file_filename = 'kernel.dat'; % file we'll use to save out kernel if precomputing

%observation_name = sprintf('%s_%s',params.observation_base,date()); % Now using the date the job started on all 
observation_path = fullfile(params.output_path_start_in,params.observation_base); % The problem with this: we're not creating the directory if it doesn't exist. That will have to change.


if params.is_L_data

    N_col = 30;
    N_row = round(N_col/num_cols*num_rows);
    length_window = round(num_bands/7)*2+1;

    [spec_pnt,spectrogram] = prelearn_pnt_params(crism_iof,N_row,N_col,length_window,params.spec_range);
    figure;plot(wa(floor(num_cols/2),:)/1000,spec_pnt);
    if params.save_spectrogram_flag
        w = [1:length_window]/length_window;
        figure;imagesc(wa(floor(num_cols/2),:)/1000,w(round(length_window/10)+1:end),log10(spectrogram(round(length_window/10)+1:end,:)));
        colorbar;
        xlabel('Wavelength (\mum)');ylabel('Frequency');
        title('Spectrogram (in log10)');
        spectrogram_savefilename = sprintf('%s_spectrogram.png',observation_path);
        saveas(gcf,spectrogram_savefilename);
    end
	
	wa_save = wa(round(num_cols/2),:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%    spatial setup     %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [c,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_positions,pix_normals,output_mask] = crism_spatial_setup(params.geoflags,params.filenames,ddr_iof,params.pix_spacing,params.pix_value,num_rows,num_cols,num_bands,params.size_kernel);
    output_mask = repmat(reshape(output_mask,[c_num_rows,c_num_cols,1]),[1,1,num_bands]);

    if params.start_from_scene_flag

        if any( size(c) ~= size(start_from_scene) )
            error('Scene to start from has different size than c ( (%i,%i,%i) vs (%i,%i,%i) )\n',...
                size(start_from_scene,1), size(start_from_scene,2),size(start_from_scene,3), ...
                size(c,1), size(c,2),size(c,3));
        end

        c = start_from_scene;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Kernel Precomputation  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %if geoflags.offnadir && geoflags.precompute
    if params.geoflags.precompute
        tic_kernel = tic;
        display('Computing kernels...');

        %calculate_spatial_kernel(kernel_file_filename,geoflags,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,pix_spacing, mro_position,pix_positions,pix_normals)
        huge_kernel = calculate_spatial_kernel_anygeo(params.geoflags,0,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
        display('Computing kernels...done');
        toc(tic_kernel)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%compute sensitivity factor%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic_sensitivity = tic;
    display(sprintf('\nComputing sensitivity image...'));

    mu = ones(num_rows, num_cols, num_bands); %cube to hold guess at measured scene

    %spectrally backproject, then spatially backproject 
    params.geoflags.forward = false;
    H_temp = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu,sb,wa,num_rows,num_cols,num_bands,params.size_kernel);
    if params.geoflags.precompute
        H = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, H_temp,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
    else
        H = calculate_spatial_kernel_anygeo(params.geoflags,H_temp,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
    end
    clearvars H_temp;

    display('Computing sensitivity image...done');
    toc(tic_sensitivity)

    Cost_function = zeros(1,params.num_iter+1);
    Data_fit = zeros(1,params.num_iter+1);
    d_Penalty_term = zeros(1,params.num_iter+1);
    o_Penalty_term = zeros(1,params.num_iter+1);
    running_time = zeros(1,params.num_iter+1);
    
    for k = 1:params.num_iter
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%spatial TF component%%
        %%%%%%%%%%%%%%%%%%%%%%%%

        tic_iter=tic;
        display(sprintf('\nStarting iteration %d...',k));
        fprintf('The time is %s',info_string(params.observation_base,version, algo_type));

        params.geoflags.forward = true;
        if params.geoflags.precompute
            mu = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, c,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
        else
            mu = calculate_spatial_kernel_anygeo(params.geoflags,c,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%spectral TF component%%
        %%%%%%%%%%%%%%%%%%%%%%%%%

        mu = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu,sb,wa,num_rows,num_cols,num_bands,params.size_kernel);

        previous_forward = mu; %save this to compute cost function;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%begin backprojecting process%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %find the error in the estimate
        ratio = crism_iof./ mu; %we want this ratio to approach 1, and we will be using this ratio to backproject

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%backproject spectral component%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        params.geoflags.forward = false;
        factor_first = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,ratio,sb,wa,num_rows,num_cols,num_bands,params.size_kernel);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%backproject spatial component%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if params.geoflags.precompute
            factor = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, factor_first,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
        else
            factor = calculate_spatial_kernel_anygeo(params.geoflags,factor_first,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
        end
        
        % penalty
        [temp, temp2, c] = compute_logcosh_penalty_H(H, spec_pnt, factor, c, params.beta_spat, params.delta_spat, params.beta_spec, params.delta_spec, params.spectral_span, params.max_iter, params.spatial_penalty_flag, params.spectral_penalty_flag, params.use_mex_versions);

        o_Penalty_term(k+1) = temp;
        d_Penalty_term(k+1) = temp2;
        Data_fit(k) = -sum(sum(sum(crism_iof.*log(previous_forward) - previous_forward)));
        Cost_function(k) = Data_fit(k) + o_Penalty_term(k);

        % Save the scene as of this iter, if desired
        if any(params.saving_iters==k)

            if params.save_c_on_iters_flag

                % NOTE: C is not actually map projected right. It is in a local
                % cartesian coordinate system.
                c_savefilename = sprintf('%s_c_iter_%i',observation_path,k);
                [c_flipped, localCartesianMapProjInfo] = crism_map_project(c, 'None', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing);
                write_envi_data(c_flipped,c_savefilename,'bsq' ...
                    , 'desc', desc_text ...
                    , 'wavelengths', wa_save ...
                    , 'default_bands', default_bands ...
                    , 'mapProjectionInfo', localCartesianMapProjInfo ...
                );
            end
            
            if params.save_proj_on_iters_flag 
                % Saving the map-projected version of C
                
                % First, limit the region we save to the projected file based on
                % no-boundaries mask
                c_save = nan(c_num_rows,c_num_cols,num_bands);
                c_save(output_mask) = c(output_mask);

                [projected_scene, mapProjectedMapProjInfo] = crism_map_project(c_save, 'Mars_Equirect', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing_output);

                proj_savefilename = sprintf('%s_proj_iter_%i',observation_path,k);
                write_envi_data(projected_scene,proj_savefilename,'bsq' ...
                    , 'desc', desc_text ...
                    , 'wavelengths', wa_save ...
                    , 'default_bands', default_bands ...
                    , 'mapProjectionInfo', mapProjectedMapProjInfo ...
                );
            end

            if params.save_mu_on_iters_flag
                mu_savefilename = sprintf('%s_mu_iter_%i',observation_path,k);
                write_envi_data(mu,mu_savefilename,'bsq' ...
                    , 'desc', desc_text ...
                    , 'wavelengths', wa_save ...
                    , 'default_bands', default_bands ...
                );
            end

        end

        % Switchover control: runs given code once at a certain iteration
        if params.switchover_flag && params.switchover_at_iter==k
            eval(params.switchover_code);
        end

        display(sprintf('\nIteration %d...done',k));
        running_time(k) = toc(tic_iter);
        fprintf('Elapsed time is %f seconds.\n\n',running_time(k));

    end

    fprintf('The time is %s',info_string(params.observation_base,version, algo_type));

    params.geoflags.forward = true;
    if params.geoflags.precompute
        mu = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, c,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
    else
        mu = calculate_spatial_kernel_anygeo(params.geoflags,c,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
    end
    mu = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu,sb,wa,num_rows,num_cols,num_bands,params.size_kernel);

    final_forward=mu;
    Data_fit(params.num_iter+1) = -sum(sum(sum(crism_iof.*log(final_forward) - final_forward)));
    Cost_function(params.num_iter+1) = Data_fit(params.num_iter+1) + o_Penalty_term(params.num_iter+1);

    % Save penalty info
    save(sprintf('%s_penalty.mat',observation_path),'Cost_function','Data_fit','d_Penalty_term','o_Penalty_term','running_time');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%. Write final versions of the scene and the measurements %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if params.save_c_on_iters_flag
        % NOTE: C is not actually map projected right. It is in a local
        % cartesian coordinate system.
        c_savefilename = sprintf('%s_c_iter_%i',observation_path,k);
        [c_flipped, localCartesianMapProjInfo] = crism_map_project(c, 'None', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing);
        write_envi_data(c_flipped,c_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
            , 'mapProjectionInfo', localCartesianMapProjInfo ...
        );
    end
    
    if params.save_proj_on_iters_flag
        % Saving the map-projected version of C
        
        % First, limit the region we save to the projected file based on
        % no-boundaries mask
        c_save = nan(c_num_rows,c_num_cols,num_bands);
        c_save(output_mask) = c(output_mask);
        [projected_scene, mapProjectedMapProjInfo] = crism_map_project(c_save, 'Mars_Equirect', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing_output);
        
        proj_savefilename = sprintf('%s_proj_iter_%i',observation_path,k);
        write_envi_data(projected_scene,proj_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
            , 'mapProjectionInfo', mapProjectedMapProjInfo ...
        );
    end
    
    if params.save_mu_on_iters_flag
        mu_savefilename = sprintf('%s_mu_iter_%i',observation_path,k);
        write_envi_data(mu,mu_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
        );
    end
    
    if params.save_input_flag
        ssa_savefilename = sprintf('%s_ssa',observation_path);
        write_envi_data(crism_iof,ssa_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
        );
		
        ddr_savefilename = sprintf('%s_DDR',observation_path);
        write_envi_data(ddr_iof,ddr_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'band_names', ddr_band_names ...
        );
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    format;

    fprintf('Finished.\n');
    fprintf('The time is %s',info_string(params.observation_base,version, algo_type));

    display(sprintf('\nTotal run time: '));
    toc(tic_total);


    %plot the convergence curve
    if params.spatial_penalty_flag || params.spectral_penalty_flag
        figure;
        hold on;
        semilogy(Cost_function(3:params.num_iter+1), 'r');
        semilogy(Data_fit(3:params.num_iter+1), 'b');
        legend('Cost function', 'Data fit negative log likelihood');
        title('Convergence curve with penalty applied');
        hold off;
        figure(2);
        hold on;
        semilogy(d_Penalty_term(3:params.num_iter+1), 'r');
        semilogy(o_Penalty_term(3:params.num_iter+1), 'b');
        legend('Decomposed penalty', 'Original penalty');
        hold off;
    else
        figure;
        hold on;
        semilogy(Cost_function(3:params.num_iter+1), 'r');
        legend('Cost function negative log likelihood');
        title('Convergence curve NO penalty applied');
        hold off;
    end

else % s data ------------------------------------------------------------------------------------------------------------------------------------

    display('Starting S data separate computation...');
    N_col = 30;
    N_row = round(N_col/num_cols*num_rows);
    length_window = round(num_bands/7)*2+1;
    params.spec_range = 2;
    [spec_pnt,spectrogram] = prelearn_pnt_params(crism_iof,N_row,N_col,length_window,params.spec_range);
    figure;plot(wa(floor(num_cols/2),:)/1000,spec_pnt);
    if params.save_spectrogram_flag
    w = [1:length_window]/length_window;
    figure;imagesc(wa(floor(num_cols/2),:)/1000,w(round(length_window/10)+1:end),log10(spectrogram(round(length_window/10)+1:end,:)));
    colorbar;
    xlabel('Wavelength (\mum)');ylabel('Frequency');
    title('Spectrogram (in log10)');
    spectrogram_savefilename = sprintf('%s_spectrogram.png',observation_path);
    saveas(gcf,spectrogram_savefilename);
    end

    % cut the data into two parts

    pre_end = params.pre_end - params.bands_bounds(1) + 1;
    post_start = params.post_start - params.bands_bounds(1) + 1;
    crism_iof_pre = crism_iof(:,:,1:pre_end);
    crism_iof_post = crism_iof(:,:,post_start:end);
    [~,~,num_pre] = size(crism_iof_pre);
    [~,~,num_post] = size(crism_iof_post);
    sb_pre = sb(:,:,1:pre_end);
    wa_pre = wa(:,1:pre_end);
    sb_post = sb(:,:,post_start:end);
    wa_post = wa(:,post_start:end);
    spec_pnt_pre = spec_pnt(1:pre_end);
    spec_pnt_post = spec_pnt(post_start:end);
    wa_save = [wa_pre(round(num_cols/2),:),wa_post(round(num_cols/2),:)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%    spatial setup     %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [c_pre,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_positions,pix_normals,output_mask] = crism_spatial_setup(params.geoflags,params.filenames,ddr_iof,params.pix_spacing,params.pix_value,num_rows,num_cols,num_pre,params.size_kernel);
    [c_post,~,~,~,~,~,~,~,~,~] = crism_spatial_setup(params.geoflags,params.filenames,ddr_iof,params.pix_spacing,params.pix_value,num_rows,num_cols,num_post,params.size_kernel);
    output_mask = repmat(reshape(output_mask,[c_num_rows,c_num_cols,1]),[1,1,num_post+num_pre]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Kernel Precomputation  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %if geoflags.offnadir && geoflags.precompute
    if params.geoflags.precompute
        tic_kernel = tic;
        display('Computing kernels...');

        %calculate_spatial_kernel(kernel_file_filename,geoflags,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,pix_spacing, mro_position,pix_positions,pix_normals)
        huge_kernel = calculate_spatial_kernel_anygeo(params.geoflags,0,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals);
        display('Computing kernels...done');
        toc(tic_kernel)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%compute sensitivity factor%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic_sensitivity = tic;
    display(sprintf('\nComputing sensitivity image...'));

    mu_pre = ones(num_rows, num_cols, num_pre); %cube to hold guess at measured scene
    mu_post = ones(num_rows, num_cols, num_post);
    
    %spectrally backproject, then spatially backproject 
    params.geoflags.forward = false;
    H_temp_pre = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu_pre,sb_pre,wa_pre,num_rows,num_cols,num_pre,params.size_kernel);
    H_temp_post = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu_post,sb_post,wa_post,num_rows,num_cols,num_post,params.size_kernel);
    if params.geoflags.precompute
        H_pre = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward,H_temp_pre,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
        H_post = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward,H_temp_post,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
    else
        H_pre = calculate_spatial_kernel_anygeo(params.geoflags,H_temp_pre,ddr_iof,ddr_lat,ddr_lon,sb_pre,num_rows,num_cols,num_pre,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
        H_post = calculate_spatial_kernel_anygeo(params.geoflags,H_temp_post,ddr_iof,ddr_lat,ddr_lon,sb_post,num_rows,num_cols,num_post,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
    end

    clearvars H_temp;

    display('Computing sensitivity image...done');
    toc(tic_sensitivity)

    Cost_function = zeros(1,params.num_iter+1);
    Data_fit = zeros(1,params.num_iter+1);
    d_Penalty_term = zeros(1,params.num_iter+1);
    o_Penalty_term = zeros(1,params.num_iter+1);
    running_time = zeros(1,params.num_iter+1);
    
    for k = 1:params.num_iter
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%spatial TF component%%
        %%%%%%%%%%%%%%%%%%%%%%%%

        tic_iter=tic;
        display(sprintf('\nStarting iteration %d...',k));
        fprintf('The time is %s',info_string(params.observation_base,version, algo_type));

        params.geoflags.forward = true;
        
        if params.geoflags.precompute
            mu_pre = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, c_pre,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
            mu_post = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, c_post,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
        else
            mu_pre = calculate_spatial_kernel_anygeo(params.geoflags,c_pre,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_pre,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
            mu_post = calculate_spatial_kernel_anygeo(params.geoflags,c_post,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_post,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%spectral TF component%%
        %%%%%%%%%%%%%%%%%%%%%%%%%

        mu_pre = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu_pre,sb,wa,num_rows,num_cols,num_pre,params.size_kernel);
        mu_post = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu_post,sb,wa,num_rows,num_cols,num_post,params.size_kernel);

        previous_forward_pre = mu_pre; %save this to compute cost function;
        previous_forward_post = mu_post; %save this to compute cost function;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%begin backprojecting process%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %find the error in the estimate
        ratio_pre = crism_iof_pre./ mu_pre; %we want this ratio to approach 1, and we will be using this ratio to backproject
        ratio_post = crism_iof_post./ mu_post; %we want this ratio to approach 1, and we will be using this ratio to backproject
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%backproject spectral component%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        params.geoflags.forward = false;
        factor_first_pre = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,ratio_pre,sb_pre,wa_pre,num_rows,num_cols,num_pre,params.size_kernel);
        factor_first_post = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,ratio_post,sb_post,wa_post,num_rows,num_cols,num_post,params.size_kernel);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%backproject spatial component%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if params.geoflags.precompute
            factor_pre = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, factor_first_pre,num_rows,num_cols,num_pre,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
            factor_post = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, factor_first_post,num_rows,num_cols,num_post,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
        else
            factor_pre = calculate_spatial_kernel_anygeo(params.geoflags,factor_first_pre,ddr_iof,ddr_lat,ddr_lon,sb_pre,num_rows,num_cols,num_pre,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
            factor_post = calculate_spatial_kernel_anygeo(params.geoflags,factor_first_post,ddr_iof,ddr_lat,ddr_lon,sb_post,num_rows,num_cols,num_post,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
        end

        % penalty
        [temp_pre, temp2_pre, c_pre] = compute_logcosh_penalty_H(H_pre, spec_pnt_pre, factor_pre, c_pre, params.beta_spat, params.delta_spat, params.beta_spec, params.delta_spec, params.spectral_span, params.max_iter, params.spatial_penalty_flag, params.spectral_penalty_flag, params.use_mex_versions);
        [temp_post, temp2_post, c_post] = compute_logcosh_penalty_H(H_post, spec_pnt_post, factor_post, c_post, params.beta_spat, params.delta_spat, params.beta_spec, params.delta_spec, params.spectral_span, params.max_iter, params.spatial_penalty_flag, params.spectral_penalty_flag, params.use_mex_versions);
        
        o_Penalty_term(k) = temp_pre + temp_post;
        d_Penalty_term(k) = temp2_pre + temp2_post;
        Data_fit(k) = -sum(sum(sum(crism_iof_pre.*log(previous_forward_pre) - previous_forward_pre))) -sum(sum(sum(crism_iof_post.*log(previous_forward_post) - previous_forward_post)));
        Cost_function(k) = Data_fit(k) + o_Penalty_term(k);
        
        % TODO: is this really still something we should do?
%         if k > start_iter + 1
%             if Cost_function(k) > Cost_function(k-1) 
%                break; 
%             end
%         end
        
        % Save the scene as of this iter, if desired
        if any(params.saving_iters==k)

            if params.save_c_on_iters_flag
                % NOTE: C is not actually map projected right. It is in a local
                % cartesian coordinate system.
                c_savefilename = sprintf('%s_c_iter_%i',observation_path,k);
                c_save = cat(3,c_pre,c_post);
                [c_flipped, localCartesianMapProjInfo] = crism_map_project(c_save, 'None', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing);
                write_envi_data(c_flipped,c_savefilename,'bsq' ...
                    , 'desc', desc_text ...
                    , 'wavelengths', wa_save ...
                    , 'default_bands', default_bands ...
                    , 'mapProjectionInfo', localCartesianMapProjInfo ...
                );
            end
            
            if params.save_proj_on_iters_flag 
                % Saving the map-projected version of C

                % First, limit the region we save to the projected file based on
                % no-boundaries mask
                c_save = nan(c_num_rows,c_num_cols,num_pre+num_post);
                tmp = cat(3,c_pre,c_post);
                c_save(output_mask) = tmp(output_mask);

                [projected_scene, mapProjectedMapProjInfo] = crism_map_project(c_save, 'Mars_Equirect', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing_output);

                proj_savefilename = sprintf('%s_proj_iter_%i',observation_path,k);
                write_envi_data(projected_scene,proj_savefilename,'bsq' ...
                    , 'desc', desc_text ...
                    , 'wavelengths', wa_save ...
                    , 'default_bands', default_bands ...
                    , 'mapProjectionInfo', mapProjectedMapProjInfo ...
                );
            end

            if params.save_mu_on_iters_flag
                mu_savefilename = sprintf('%s_mu_iter_%i',observation_path,k);
                mu_save = cat(3,mu_pre,mu_post);
                write_envi_data(mu_save,mu_savefilename,'bsq' ...
                    , 'desc', desc_text ...
                    , 'wavelengths', wa_save ...
                    , 'default_bands', default_bands ...
                );
            end
            
        end

        % Switchover control: runs given code once at a certain iteration
        if params.switchover_flag && params.switchover_at_iter==k
            eval(params.switchover_code);
        end

        display(sprintf('\nIteration %d...done',k));
        running_time(k) =  toc(tic_iter);
        fprintf('Elapsed time is %f seconds.\n\n',running_time(k));

    end

    fprintf('The time is %s',info_string(params.observation_base,version, algo_type));

    params.geoflags.forward = true;
    if params.geoflags.precompute
        mu_pre = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, c_pre,num_rows,num_cols,num_pre,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
        mu_post = general_proj_spat_precomputed(huge_kernel, params.geoflags.forward, c_post,num_rows,num_cols,num_post,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
    else
        mu_pre = calculate_spatial_kernel_anygeo(params.geoflags,c_pre,ddr_iof,ddr_lat,ddr_lon,sb_pre,num_rows,num_cols,num_pre,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
        mu_post = calculate_spatial_kernel_anygeo(params.geoflags,c_post,ddr_iof,ddr_lat,ddr_lon,sb_post,num_rows,num_cols,num_post,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
    end
    mu_pre = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu_pre,sb_pre,wa_pre,num_rows,num_cols,num_pre,params.size_kernel);
    mu_post = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu_post,sb_post,wa_post,num_rows,num_cols,num_post,params.size_kernel);

    final_forward_pre=mu_pre;
    final_forward_post=mu_post;
    Data_fit(params.num_iter+1) = -sum(sum(sum(crism_iof_pre.*log(final_forward_pre) - final_forward_pre))) -sum(sum(sum(crism_iof_post.*log(final_forward_post) - final_forward_post)));
    Cost_function(params.num_iter+1) = Data_fit(params.num_iter+1) + o_Penalty_term(params.num_iter+1);

    % Save penalty info
    save(sprintf('%s_penalty.mat',observation_path),'Cost_function','Data_fit','d_Penalty_term','o_Penalty_term','running_time');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%. Write final versions of the scene and the measurements %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if params.save_c_on_iters_flag
        % NOTE: C is not actually map projected right. It is in a local
        % cartesian coordinate system.
        c_savefilename = sprintf('%s_c_iter_%i',observation_path,k);
        c_save = cat(3,c_pre,c_post);
        [c_flipped, localCartesianMapProjInfo] = crism_map_project(c_save, 'None', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing);
        write_envi_data(c_flipped,c_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
            , 'mapProjectionInfo', localCartesianMapProjInfo ...
        );
    end

    if params.save_proj_on_iters_flag 
        % Saving the map-projected version of C

        % First, limit the region we save to the projected file based on
        % no-boundaries mask
        c_save = nan(c_num_rows,c_num_cols,num_pre+num_post);
        tmp = cat(3,c_pre,c_post);
        c_save(output_mask) = tmp(output_mask);

        [projected_scene, mapProjectedMapProjInfo] = crism_map_project(c_save, 'Mars_Equirect', ddr_iof, ddr_lat, ddr_lon, params.pix_spacing_output);

        proj_savefilename = sprintf('%s_proj_iter_%i',observation_path,k);
        write_envi_data(projected_scene,proj_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
            , 'mapProjectionInfo', mapProjectedMapProjInfo ...
        );
    end

    if params.save_mu_on_iters_flag
        mu_savefilename = sprintf('%s_mu_iter_%i',observation_path,k);
        mu_save = cat(3,mu_pre,mu_post);
        write_envi_data(mu_save,mu_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
            );
    end
    
    % Saves that only occur at the end
    if params.save_input_flag
        ssa_save = cat(3,crism_iof_pre,crism_iof_post);
        ssa_savefilename = sprintf('%s_ssa',observation_path);
        write_envi_data(ssa_save,ssa_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'wavelengths', wa_save ...
            , 'default_bands', default_bands ...
        );

        ddr_savefilename = sprintf('%s_DDR',observation_path);
        write_envi_data(ddr_iof,ddr_savefilename,'bsq' ...
            , 'desc', desc_text ...
            , 'band_names', ddr_band_names ...
        );
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    format;

    fprintf('Finished.\n');
    fprintf('The time is %s',info_string(params.observation_base,version, algo_type));

    display(sprintf('\nTotal run time: '));
    toc(tic_total);


    %plot the convergence curve
    if params.spatial_penalty_flag || params.spectral_penalty_flag
        figure;
        hold on;
        semilogy(Cost_function(3:params.num_iter+1), 'r');
        semilogy(Data_fit(3:params.num_iter+1), 'b');
        legend('Cost function', 'Data fit negative log likelihood');
        title('Convergence curve with penalty applied');
        hold off;
        figure;
        hold on;
        semilogy(d_Penalty_term(3:params.num_iter+1), 'r');
        semilogy(o_Penalty_term(3:params.num_iter+1), 'b');
        legend('Decomposed penalty', 'Original penalty');
        hold off;
    else
        figure;
        hold on;
        semilogy(Cost_function(3:params.num_iter+1), 'r');
        legend('Cost function negative log likelihood');
        title('Convergence curve NO penalty applied');
        hold off;
    end

end
end

function [info] = info_string(obs_name,version,algo)
% INFO_STRING returns string that identifies this run in the printed logs
    info = sprintf('%s (scene: %s; version: %s %s)\n',datestr(now),obs_name,algo,version);
end
