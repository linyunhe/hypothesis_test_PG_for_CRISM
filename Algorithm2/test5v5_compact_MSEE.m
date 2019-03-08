% Now encapsulated in a function to
%   1) Make it runnable on CHPC cluster,
%   2) Prevent polution of the namespace
%
% It can have no arguments, so it can still be executed by pressing "run."
%
% ARGUMENTS: Takes pairs of arguments of the form
% 'ParameterName','ParameterValue'. The name should be that of some field
% of params found at the top of the function, and the value should be of
% whatever matlab type is appropriate. If argume nts are given, they are
% used instead of whatever values are written into the script; the script
% values can be considered the default values.
function [Cost_function] = test5v5_compact_MSEE(varargin)

format;

tic_total = tic;

%%%%%%%%%%%%%
%%edit this%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.spatial_penalty_flag = 0; % 1 for using penalty, 0 for no penalty, except 20 on no penalty 20 on with penalty
params.spectral_penalty_flag = 0; % 1 for using penalty, 0 for no penalty, except 20 on with no penalty 20 on with penalty
params.beta_spat = 0.02;
params.delta_spat = 4.0;
params.beta_spec = 0.2;
params.delta_spec = 0.9;
params.max_iter = 3; %maximum iterations be used in Restricted Region Newton's method
params.spectral_span = 5; %number of spectral neighbours considered for spectral penalty

params.use_mex_versions = true; % Use MEX C++ translations when possible
params.mex_num_threads = -1; % Number of threads to request. -1 requests the maximum number.

global mex_num_threads; mex_num_threads = params.mex_num_threads;
% Looking for num_rows, num_cols, num_bands? No longer necessary, since we
% can read these from data header files

params.is_L_data = 1; % true for L, false for S. For choosing spectral transfer method

params.size_kernel = 11; %must be odd
params.pix_spacing = 12.0; %size of pixel diameter in desired image
params.pix_value = 0.4; %default intial c pixel value. Must be greater than 0 or less than 1.
params.num_iter = 250; %number of times to iterate through the algorithm.
params.geoflags.precompute = 0;
params.g1_flag = true; % if true, use 1 spectral gaussian even for L data, otherwise use 3 assym. gaus. in that case

%------------Spatial TF parameters-------------------------
params.geoflags.offnadir = 0; % If true: use flat, off-nadir geo. If false: use nadir geo.
%The following are only applicable if using new geometry
params.geoflags.flat = 1; %Approximating a flat surface (still accounting for viewing angle)
%-----------------------------------------------------------
% display('1');
% Switchover feature. If true, the given fragment of code is executed when the given iteration finishes
params.switchover_flag = 1;
params.switchover_at_iter = 5; % num of the iteration after which we apply the switchover
params.switchover_code = 'params.spatial_penalty_flag=1;params.spectral_penalty_flag=1;';

params.spectral_offset = -0.5294; % in units nm

params.save_c_on_iters_flag = 1; % whether to save the scene on the iters listed below
params.save_mu_on_iters_flag = 0; % whether to save sensor space on the iters listed below
params.saving_iters = [];
% display('2');
% input files
params.crism_iof_filename = 'E:\simulations_MSEE_MLM\Simulations2.bsq'; %derived SSA alreadymedian filtered
params.ddr_iof_filename = 'E:\06_06_2016_pangboche\paper\Hypothesis_simulations\ATO00037D74_01_DE169L_DDR1_envi.img'; %gives latitude, logitude, emission angle, etc
params.sb_filename = 'F:/data/sb_and_wa/cdr490947778566_sb0000000l_3_envi.img'; %parameters for spectral TF
params.wa_filename = 'F:/data/sb_and_wa/cdr490947778566_wa0000000l_3_envi.img'; %bandpass centers

params.filenames.radius_filename = '';
params.filenames.topo_filename = '';
params.filenames.spice_filename = '';
params.use_mex_versions= 1;
params.start_from_scene_flag = 0;
params.start_from_scene_filename = ''; % an ENVI-format scene that we use as initial guess

% Paths to prepend to input and output file paths respectively ("." is the
% current directory). If empty string, assumes an acceptable full path has
% been given.
% display('3');
params.input_path_start_in = '';
params.output_path_start_in = '';

% input file subsetting
% (all indexed from 1)
% if the *_use_all flag is false, uses the *_subset_bounds as the min and max extent
% only datasets with 640 columns are supported (and matching # of rows & bands in SSA & DDR)
params.rows_use_all = false;
params.rows_subset_bounds = [1,170];
params.cols_use_all = true;
params.cols_subset_bounds = [250,350];
params.bands_bounds = [193,197]; % the bands that corresond to the range found in the SSA

% output file naming
params.observation_base = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use function args to replace some params %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display('4');
params = setParamsFromArgs( varargin, params );

% build a description we can use in header files
desc_text = generateParamListDesc(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data loading from the files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display('5');    
% read data files
crism_iof = multibandread_fromENVIHeader(fullfile(params.input_path_start_in,params.crism_iof_filename)); %measured radiance
ddr_iof = multibandread_fromENVIHeader(fullfile(params.input_path_start_in,params.ddr_iof_filename)); %gives latitude, logitude, emission angle, etc
sb = multibandread_fromENVIHeader(fullfile(params.input_path_start_in,params.sb_filename)); %parameters for spectral TF
wa = multibandread_fromENVIHeader(fullfile(params.input_path_start_in,params.wa_filename)); %bandpass centers
params.filenames.radius_filename = fullfile(params.input_path_start_in,params.filenames.radius_filename);
params.filenames.topo_filename = fullfile(params.input_path_start_in,params.filenames.topo_filename);
params.filenames.spice_filename = fullfile(params.input_path_start_in,params.filenames.spice_filename);
% crism_iof = repmat(crism_iof,[178 ,1,1]);
% ddr_iof = ddr_iof(1:170,:,:);
% crism_iof = crism_iof(1:170,:,:);
% display('6');

if params.start_from_scene_flag
    % An ENVI-format scene that we use as initial guess.
    % Currently must be the same size as the scene crism_setup would create
    % for us
    start_from_scene = multibandread_fromENVIHeader(fullfile(params.input_path_start_in,params.start_from_scene_filename));
	start_from_scene = flipdim(flipdim(start_from_scene,1),2); % flip over both spatial axes to restore to expected orientation
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
elseif size(crism_iof,2) ~= 640
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
    error('IOF/SSA and DDR must have same number of rows (currently, %i and %i)',size(crism_iof,1),size(ddr_iof,1));
elseif size(sb,3)~=size(wa,3)
    error('SB and WA must have same number of bands (currently, %i and %i)',size(sb,3),size(wa,3));
end

if params.is_L_data
    col_furthest_bounds = [32,631];
    params.spec_range = 10;
else % S data
    col_furthest_bounds = [26,626];
    params.spec_range = 2;
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
    col_min = params.rows_subset_bounds(1);
    col_max = params.rows_subset_bounds(2);
end

band_min = params.bands_bounds(1);
band_max = params.bands_bounds(2);

% make sure no columns requested that are out-of-bounds
if (col_min < col_furthest_bounds(1)) || (col_furthest_bounds(2) < col_max)
    error('Invalid columns selected [%i, %i]. Must be within (inclusive) range [%i, %i].', ...
        col_min,col_max, ...
        col_furthest_bounds(1), col_furthest_bounds(2));
end

% one last consistency check: is the correct number of bands selected?
% band_min = 1;
% band_max = 246;

crism_iof(crism_iof < 0) = 0;
if band_max - band_min + 1 ~= size(crism_iof,3)
    error('The selected band range [%i, %i] does not have the correct number of bands (%i; currently: %i)',band_min,band_max,size(crism_iof,3),band_max-band_min+1);
end

% currently, no support for subsetting bands within SSA
% row_min = 100;
% row_max = 130;
% col_min = 300;
% col_max = 360;
crism_iof = crism_iof(row_min:row_max,col_min:col_max,:); 
ddr_iof = ddr_iof(row_min:row_max,col_min:col_max,:);
sb = sb(:,col_min:col_max,band_min:band_max);
wa = wa(:,col_min:col_max,band_min:band_max);

% capture dimension sizes in variables that the code uses
[num_rows, num_cols, num_bands] = size(crism_iof);

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

fprintf('Starting run on scene with dimensions %i rows, %i cols, %i bands\n',num_rows,num_cols,num_bands);

kernel_file_filename = 'kernel.dat'; % file we'll use to save out kernel if precomputing

%observation_name = sprintf('%s_%s',params.observation_base,date()); % Now using the date the job started on all 
observation_path = fullfile(params.output_path_start_in,params.observation_base); % The problem with this: we're not creating the directory if it doesn't exist. That will have to change.
spec_pnt = 0;
if  params.switchover_flag || params.spectral_penalty_flag
    N_col = 30;
    N_row = round(N_col/num_cols*num_rows);
    length_window = round(num_bands/7)*2+1;

%     [spec_pnt,spectrogram] = prelearn_pnt_params(crism_iof,N_row,N_col,length_window,params.spec_range);
    spec_pnt = 5*ones(1,num_bands);
%     figure;plot(wa(floor(num_cols/2),:)/1000,spec_pnt);
%         w = [1:length_window]/length_window;
%         figure;imagesc(wa(floor(num_cols/2),:)/1000,w(round(length_window/10)+1:end),log10(spectrogram(round(length_window/10)+1:end,:)));
%         colorbar;
%         xlabel('Wavelength (\mum)');ylabel('Frequency');
%         title('Spectrogram (in log10)');
%         spectrogram_savefilename = sprintf('%s_spectrogram.png',observation_path);
%         saveas(gcf,spectrogram_savefilename);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    spatial setup     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[c,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_positions,pix_normals] = crism_spatial_setup(params.geoflags,params.filenames,ddr_iof,params.pix_spacing,params.pix_value,num_rows,num_cols,num_bands,params.size_kernel);

if params.start_from_scene_flag
    
    if any( size(c) ~= size(start_from_scene) )
        error('Scene to start from has different size than c ( (%i,%i,%i) vs (%i,%i,%i) )\n',...
            size(start_from_scene,1), size(start_from_scene,2),size(start_from_scene,3), ...
            size(c,1), size(c,2),size(c,3));
    end
    
    c = start_from_scene;
end

% weight = 1/(11*11*11+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kernel Precomputation  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if geoflags.newgeo && geoflags.precompute
if params.geoflags.precompute
    tic_kernel = tic;
    display('Computing kernels...');
    
    %calculate_spatial_kernel(kernel_file_filename,geoflags,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,pix_spacing, mro_position,pix_positions,pix_normals)
    f = calculate_spatial_kernel_anygeo(kernel_file_filename,params.geoflags,0,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals);
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
% Is_H = true;
H_temp = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu,sb,wa,num_rows,num_cols,num_bands,params.size_kernel);
if params.geoflags.precompute
    H = general_proj_spat_precomputed(kernel_file_filename, params.geoflags.forward, H_temp,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
else
    H = calculate_spatial_kernel_anygeo(params.geoflags,H_temp,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
end
% ratio = 0.9;
% weight = H*ratio;
% H = H./weight;
clearvars H_temp;
display('Computing sensitivity image...done');
toc(tic_sensitivity)

Cost_function = zeros(1,params.num_iter+1);
Data_fit = zeros(1,params.num_iter+1);
d_Penalty_term = zeros(1,params.num_iter+1);
o_Penalty_term = zeros(1,params.num_iter+1);

for k = 1:params.num_iter
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%spatial TF component%%
    %%%%%%%%%%%%%%%%%%%%%%%%

    tic_iter=tic;
    display(sprintf('\nStarting iteration %d...',k));

    params.geoflags.forward = true;
    if params.geoflags.precompute
        mu = general_proj_spat_precomputed(kernel_file_filename, params.geoflags.forward, c,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
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
    err = crism_iof - mu; %we want this ratio to approach 1, and we will be using this ratio to backproject
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%backproject spectral component%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.geoflags.forward = false;
    factor_first = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,err,sb,wa,num_rows,num_cols,num_bands,params.size_kernel);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%backproject spatial component%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if params.geoflags.precompute
        factor = general_proj_spat_precomputed(kernel_file_filename, params.geoflags.forward, factor_first,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
    else
        factor = calculate_spatial_kernel_anygeo(params.geoflags,factor_first,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
    end
    
%     if k >= 6
%           a = 1;
%     end
    %solve the convex decomposed penalty function using Newton's method
    [temp, temp2, c] = compute_logcosh_penalty_H_MSEE (H, spec_pnt,factor, c, params.beta_spat, params.delta_spat, params.beta_spec, params.delta_spec, params.spectral_span, params.max_iter, params.spatial_penalty_flag, params.spectral_penalty_flag,1);
    
    o_Penalty_term(k+1) = temp;
    d_Penalty_term(k+1) = temp2;
    Data_fit(k) = sum(sum(sum((crism_iof - previous_forward).^2)));
    Cost_function(k) = Data_fit(k) + d_Penalty_term(k);
    
    % Save the scene as of this iter, if desired
    if any(params.saving_iters==k)
        if params.save_c_on_iters_flag
            c_savefilename = sprintf('%s_c_iter_%i',observation_path,k);
            [flippedC, mapProjectionInfo] = transformSceneAndGetMapProj( c, params.pix_spacing, ddr_lat, ddr_lon );
            multibandwrite(flippedC, sprintf('%s.bsq',c_savefilename), 'bsq');
            save_envi_header_mapproject(flippedC, sprintf('%s.hdr',c_savefilename), 'bsq', desc_text, wa(floor(num_cols/2),:), mapProjectionInfo);
        end      
        if params.save_mu_on_iters_flag
            mu_savefilename = sprintf('%s_mu_iter_%i',observation_path,k);
            multibandwrite(mu, sprintf('%s.bsq',mu_savefilename), 'bsq');
            save_envi_header(mu, sprintf('%s.hdr',mu_savefilename), 'bsq', desc_text, wa(floor(num_cols/2),:));
        end
    end
	
	% Switchover control: runs given code once at a certain iteration
	if params.switchover_flag && params.switchover_at_iter==k
		eval(params.switchover_code);
	end
    
    display(sprintf('\nIteration %d...done',k));
    toc(tic_iter);
     
end

% Final forward
% if geoflags.newgeo
%     geoflags.forward = true;
%     if geoflags.precompute
%         mu = general_proj_spat_precomputed(kernel_file_filename, geoflags.forward, c,num_rows,num_cols,num_bands,size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
%     else
%         error('Currently can''t do new geometry without precompute');
%     end
% else
%     mu = forward_spat(c,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,pix_spacing);
% end
params.geoflags.forward = true;
if params.geoflags.precompute
    mu = general_proj_spat_precomputed(kernel_file_filename, params.geoflags.forward, c,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols);
else
    mu = calculate_spatial_kernel_anygeo(params.geoflags,c,ddr_iof,ddr_lat,ddr_lon,sb,num_rows,num_cols,num_bands,params.size_kernel,lat_number,lon_number,c_num_rows,c_num_cols,params.pix_spacing, mro_position,pix_positions,pix_normals,params.use_mex_versions);
end
mu = general_proj_spec(params.is_L_data,params.geoflags.forward,params.g1_flag,mu,sb,wa,num_rows,num_cols,num_bands,params.size_kernel);

final_forward=mu;
Data_fit(params.num_iter+1) = -sum(sum(sum(crism_iof - final_forward)));
Cost_function(params.num_iter+1) = Data_fit(params.num_iter+1) + d_Penalty_term(params.num_iter+1);

% Save penalty info
%save(sprintf('%s_penalty.mat',observation_path),'Cost_function','Data_fit','d_Penalty_term','o_Penalty_term');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%. Write final versions of the scene and the measurements %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_savefilename = sprintf('%s_c_iter_%i',observation_path,k);
mu_savefilename = sprintf('%s_mu_iter_%i',observation_path,k);

[flippedC, mapProjectionInfo] = transformSceneAndGetMapProj( c, params.pix_spacing, ddr_lat, ddr_lon );
multibandwrite(flippedC, sprintf('%s.bsq',c_savefilename), 'bsq');
save_envi_header_mapproject(c, sprintf('%s.hdr',c_savefilename), 'bsq', desc_text, wa(floor(num_cols/2),:), mapProjectionInfo);

multibandwrite(mu, sprintf('%s.bsq',mu_savefilename), 'bsq');
save_envi_header(mu, sprintf('%s.hdr',mu_savefilename), 'bsq', desc_text, wa(floor(num_cols/2),:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format;

display(sprintf('\nTotal run time: '));
toc(tic_total);


%plot the convergence curve
if params.spatial_penalty_flag || params.spectral_penalty_flag
    figure;
    hold on;
    plot(Cost_function(3:params.num_iter+1), 'r');
    plot(Data_fit(3:params.num_iter+1), 'b');
    legend('Cost function', 'Data fit negative log likelihood');
    title('Convergence curve with penalty applied');
    hold off;
    figure;
    hold on;
    plot(d_Penalty_term(3:params.num_iter+1), 'r');
    plot(o_Penalty_term(3:params.num_iter+1), 'b');
    legend('Decomposed penalty', 'Original penalty');
    hold off;
else
    figure;
    hold on;
    plot(Cost_function(3:params.num_iter+1), 'r');
    legend('Cost function negative log likelihood');
    title('Convergence curve NO penalty applied');
    hold off;
end

end
