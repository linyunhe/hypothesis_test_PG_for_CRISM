function gui_react ( hObject, handles )
%GUI_REACT Manipulate GUI components based on callbacks to other components
%   hObject is the object whose callback just activated, and handles is the
%   set of all CRISM GUI handles

% No reactions for events in edit boxes
if strcmp(hObject.Style,'edit')
    return;
end

% Based on object user interacted with, call functions to manipulate the
% objects that object could change the state of
switch hObject
    % Checkboxes for "all" along spatial/spectral dimensions
    case handles.ssa_row_all
        manipulate_in_response(handles.ssa_row_min,handles);
        manipulate_in_response(handles.ssa_row_max,handles);
    case handles.ssa_col_all
        manipulate_in_response(handles.ssa_col_min,handles);
        manipulate_in_response(handles.ssa_col_max,handles);
    % Input Panel browse controls
    case handles.ssa_browse
        manipulate_in_response(handles.ssa_name,handles);
    case handles.ddr_browse
        manipulate_in_response(handles.ddr_name,handles);
    case handles.sb_browse
        manipulate_in_response(handles.sb_name,handles);
    case handles.wa_browse
        manipulate_in_response(handles.wa_name,handles);
    case handles.auto_find_sb_wa_button
        % TODO: I don't like how this breaks normal protocol in this system
        auto_fill_sb_wa(handles);
    % General Geometry Panel enables
    case handles.projected_start_flag
        manipulate_in_response(handles.proj_scene_name,handles);
        manipulate_in_response(handles.proj_scene_browse,handles);
        manipulate_in_response(handles.init_ssa_guess,handles);
    case handles.l_band_button
        manipulate_in_response(handles.g1_checkbox,handles);
        manipulate_in_response(handles.pre_end,handles);
        manipulate_in_response(handles.post_start,handles);
        manipulate_in_response(handles.spec_range,handles);
    case handles.s_band_button
        manipulate_in_response(handles.g1_checkbox,handles);
        manipulate_in_response(handles.pre_end,handles);
        manipulate_in_response(handles.post_start,handles);
        manipulate_in_response(handles.spec_range,handles);
    case handles.g1_checkbox
        % pass
    % Iterations Panel enables
    case handles.no_switch_radio
        manipulate_in_response(handles.switchover_at_iter,handles);
        manipulate_in_response(handles.switch_penalty_spectral,handles);
        manipulate_in_response(handles.switch_penalty_spatial,handles);
        % could also activate the penalty settings when you change this
        manipulate_in_response(handles.beta_spec,handles);
        manipulate_in_response(handles.delta_spec,handles);
        manipulate_in_response(handles.spec_span,handles);
        manipulate_in_response(handles.spec_range,handles);
        manipulate_in_response(handles.pen_internal_iters,handles);
        manipulate_in_response(handles.beta_spat,handles);
        manipulate_in_response(handles.delta_spat,handles);
    case handles.switch_penalty_radio
        manipulate_in_response(handles.switchover_at_iter,handles);
        manipulate_in_response(handles.switch_penalty_spectral,handles);
        manipulate_in_response(handles.switch_penalty_spatial,handles);
        % could also activate the penalty settings when you turn it on
        manipulate_in_response(handles.beta_spec,handles);
        manipulate_in_response(handles.delta_spec,handles);
        manipulate_in_response(handles.spec_span,handles);
        manipulate_in_response(handles.spec_range,handles);
        manipulate_in_response(handles.pen_internal_iters,handles);
        manipulate_in_response(handles.beta_spat,handles);
        manipulate_in_response(handles.delta_spat,handles);
    case handles.switch_penalty_spectral
        manipulate_in_response(handles.beta_spec,handles);
        manipulate_in_response(handles.delta_spec,handles);
        manipulate_in_response(handles.spec_span,handles);
        manipulate_in_response(handles.pen_internal_iters,handles);
    case handles.switch_penalty_spatial
        manipulate_in_response(handles.beta_spat,handles);
        manipulate_in_response(handles.delta_spat,handles);
        manipulate_in_response(handles.pen_internal_iters,handles);
    % Save Data Panel enables
    case handles.save_mu_flag
        manipulate_in_response(handles.save_iters,handles);
    case handles.save_c_flag
        manipulate_in_response(handles.save_iters,handles);
    case handles.save_proj_flag
        manipulate_in_response(handles.save_iters,handles);
    case handles.save_input_flag
        manipulate_in_response(handles.save_iters,handles);
    case handles.save_spectrogram_flag
        manipulate_in_response(handles.save_iters,handles);
    case handles.output_dir_browse
        manipulate_in_response(handles.output_dir,handles);
    case handles.output_dir_auto_button
        manipulate_in_response(handles.output_dir,handles);
    % Decomposed Penalty Panel enables
    case handles.spec_pen_flag
        manipulate_in_response(handles.beta_spec,handles);
        manipulate_in_response(handles.delta_spec,handles);
        manipulate_in_response(handles.spec_span,handles);
        manipulate_in_response(handles.pen_internal_iters,handles);
    case handles.spat_pen_flag
        manipulate_in_response(handles.beta_spat,handles);
        manipulate_in_response(handles.delta_spat,handles);
        manipulate_in_response(handles.pen_internal_iters,handles);
    % Geometry Selection Panel enables
    case handles.nadir_geo_button
        manipulate_in_response(handles.radius_name,handles);
        manipulate_in_response(handles.topo_name,handles);
        manipulate_in_response(handles.spice_name,handles);
        manipulate_in_response(handles.radius_browse,handles);
        manipulate_in_response(handles.topo_browse,handles);
        manipulate_in_response(handles.spice_browse,handles);
    case handles.off_nadir_flat_geo_button
        manipulate_in_response(handles.radius_name,handles);
        manipulate_in_response(handles.topo_name,handles);
        manipulate_in_response(handles.spice_name,handles);
        manipulate_in_response(handles.radius_browse,handles);
        manipulate_in_response(handles.topo_browse,handles);
        manipulate_in_response(handles.spice_browse,handles);
    case handles.off_nadir_rolling_geo_button
        manipulate_in_response(handles.radius_name,handles);
        manipulate_in_response(handles.topo_name,handles);
        manipulate_in_response(handles.spice_name,handles);
        manipulate_in_response(handles.radius_browse,handles);
        manipulate_in_response(handles.topo_browse,handles);
        manipulate_in_response(handles.spice_browse,handles);
    % Non-standard file browsing
    case handles.proj_scene_browse
        manipulate_in_response(handles.proj_scene_name,handles);
    case handles.radius_browse
        manipulate_in_response(handles.radius_name,handles);
    case handles.topo_browse
        manipulate_in_response(handles.topo_name,handles);
    case handles.spice_browse
        manipulate_in_response(handles.spice_name,handles);
    % Options file save/load
    case handles.save_button
        sel_file = browse_for_file_return_path('Save Options File',true,{'*.opt','Option files';'*.*','All files'});
        if sel_file ~= 0
            save_options_file(handles,sel_file);
        end
    case handles.load_button
        sel_file = browse_for_file_return_path('Load Options File',false,{'*.opt','Option files';'*.*','All files'});
        if sel_file ~= 0
            restore_options_file(handles,sel_file);
        end
    % Start Button Control: go into routine to parse args and start
    % algorithm!
    case handles.start_button
        start_algo(handles);
    % Otherwise, give message about nothing to do here
    otherwise
        fprintf('No components'' state depends on: %s\n',hObject.Tag);
end

end

function manipulate_in_response( obj, handles )
% Chance state of an object based on the state of others (one of which has
% changed)

switch obj
    % SSA min/max boxes
    case handles.ssa_row_min
        toggle_control( ~ get(handles.ssa_row_all,'Value'), obj);
    case handles.ssa_row_max
        toggle_control( ~ get(handles.ssa_row_all,'Value'), obj);
    case handles.ssa_col_min
        toggle_control( ~ get(handles.ssa_col_all,'Value'), obj);
    case handles.ssa_col_max
        toggle_control( ~ get(handles.ssa_col_all,'Value'), obj);
    case handles.ssa_band_min
        toggle_control( ~ get(handles.ssa_band_all,'Value'), obj);
    case handles.ssa_band_max
        toggle_control( ~ get(handles.ssa_band_all,'Value'), obj);
    % Input panel file names browsing
    case handles.ssa_name
        if get(handles.ssa_browse, 'value')
            browse_for_file( obj , 'SSA File');
        end
    case handles.ddr_name
        if get(handles.ddr_browse, 'value')
            browse_for_file( obj , 'DDR File');
        end
    case handles.sb_name
        if get(handles.sb_browse, 'value')
            browse_for_file( obj , 'SB File');
        end
    case handles.wa_name
        if get(handles.wa_browse, 'value')
            browse_for_file( obj , 'WA File');
        end
    % General geometry
    case handles.proj_scene_name
        toggle_control( get(handles.projected_start_flag,'Value'), obj);
        if get(handles.proj_scene_browse, 'value')
            browse_for_file( obj , 'Projected Scene File' );
        end
    case handles.proj_scene_browse
        toggle_control( get(handles.projected_start_flag,'Value'), obj);
    case handles.init_ssa_guess
        toggle_control( ~ get(handles.projected_start_flag,'Value'), obj);
    case handles.g1_checkbox
        toggle_control( get(handles.l_band_button, 'value'), obj);    
    case handles.pre_end
        toggle_control( get(handles.s_band_button,'value'),obj);
    case handles.post_start
        toggle_control( get(handles.s_band_button, 'value'), obj);
    case handles.spec_range
        toggle_control( get(handles.l_band_button, 'value'), obj);
    % Iterations Panel
    case handles.switchover_at_iter
        % we assume that the buttons for switchover & no switchover are
        % mutually exclusive (which is currently enforced via button group)
        toggle_control( get(handles.switch_penalty_radio,'Value'), obj);
    case handles.switch_penalty_spectral
        toggle_control( get(handles.switch_penalty_radio,'Value'), obj);
    case handles.switch_penalty_spatial
        toggle_control( get(handles.switch_penalty_radio,'Value'), obj);
    % Save Data Panel
    case handles.save_iters
        % if anything is being saved each iter, allow entering iter nums
        toggle_control( get(handles.save_c_flag,'Value') || get(handles.save_mu_flag,'Value') || get(handles.save_proj_flag,'Value') || get(handles.save_input_flag,'Value') || get(handles.save_spectrogram_flag,'Value') ...
            , obj);
    case handles.output_dir
        if get(handles.output_dir_browse, 'value')
            browse_for_dir( obj, 'Output Directory' );
        elseif get(handles.output_dir_auto_button, 'value')
            % if SSA's location has been specified, copy the directory name
            % to this box
            copy_dir_from_another_input(obj,handles.ssa_name);
        end
    % Decomposed Penalty Panel
    case handles.beta_spec
        should_show = bool_from_box(handles.spec_pen_flag) || (bool_from_box(handles.switch_penalty_spectral) && bool_from_box(handles.switch_penalty_radio));
        toggle_control( should_show, obj);
    case handles.delta_spec
        should_show = bool_from_box(handles.spec_pen_flag) || (bool_from_box(handles.switch_penalty_spectral) && bool_from_box(handles.switch_penalty_radio));
        toggle_control( should_show, obj);
    case handles.beta_spat
        should_show = bool_from_box(handles.spat_pen_flag) || (bool_from_box(handles.switch_penalty_spatial) && bool_from_box(handles.switch_penalty_radio));
        toggle_control( should_show, obj);
    case handles.delta_spat
        should_show = bool_from_box(handles.spat_pen_flag) || (bool_from_box(handles.switch_penalty_spatial) && bool_from_box(handles.switch_penalty_radio));
        toggle_control( should_show, obj);
    case handles.spec_span
        should_show = bool_from_box(handles.spec_pen_flag) || (bool_from_box(handles.switch_penalty_spectral) && bool_from_box(handles.switch_penalty_radio));
        toggle_control( should_show, obj);
    case handles.pen_internal_iters
        should_show = bool_from_box(handles.spec_pen_flag) ||  bool_from_box(handles.spat_pen_flag) ...
            || ( (bool_from_box(handles.switch_penalty_spectral) || bool_from_box(handles.switch_penalty_spatial)) && bool_from_box(handles.switch_penalty_radio) ...
        );
        toggle_control( should_show, obj);
    % Geometry Selection Panel
    case handles.radius_name
        toggle_control( get(handles.off_nadir_flat_geo_button,'Value') || get(handles.off_nadir_rolling_geo_button,'Value'), obj);
        if get(handles.radius_browse, 'value')
            browse_for_file( obj , 'radius File');
        end
    case handles.topo_name
        toggle_control( get(handles.off_nadir_flat_geo_button,'Value') || get(handles.off_nadir_rolling_geo_button,'Value'), obj);
        if get(handles.topo_browse, 'value')
            browse_for_file( obj , 'topo File');
        end
    case handles.spice_name
        toggle_control( get(handles.off_nadir_flat_geo_button,'Value') || get(handles.off_nadir_rolling_geo_button,'Value'), obj);
        if get(handles.spice_browse, 'value')
            browse_for_file( obj , 'SPICE File');
        end
    case handles.radius_browse
        toggle_control( get(handles.off_nadir_flat_geo_button,'Value') || get(handles.off_nadir_rolling_geo_button,'Value'), obj);
    case handles.topo_browse
        toggle_control( get(handles.off_nadir_flat_geo_button,'Value') || get(handles.off_nadir_rolling_geo_button,'Value'), obj);
    case handles.spice_browse
        toggle_control( get(handles.off_nadir_flat_geo_button,'Value') || get(handles.off_nadir_rolling_geo_button,'Value'), obj);
    % Otherwise, give message about nothing to do here
    otherwise
        fprintf('Component has no instructions to change state: %s\n',obj.Tag);
end

end

% Enable/disable a given control based on a bool value
function toggle_control(newstate, ctrl)
% newstate is a bool; we enable if true and disable if false
    if newstate
        set(ctrl,'enable','on');
    else
        set(ctrl,'enable','off');
    end
end

% Browse for a file and enter it in given edit object
function browse_for_file(textBox,titleStr)
    global last_browsed_path;
    [filename, pathname, ~] = uigetfile({'*.*','All files'},titleStr,last_browsed_path);
    fullpath = fullfile(pathname, filename);
    if filename ~= 0 % user did not press cancel
        set(textBox,'String',fullpath);
        last_browsed_path = pathname;
    end
end

% Browse for a file an return its path for later use, returning 0 for null
function selected_file = browse_for_file_return_path(titleStr,is_saving,filetype_cells)
    global last_browsed_path;
    if is_saving
        [filename, pathname, ~] = uiputfile(filetype_cells,titleStr,last_browsed_path);
    else
        [filename, pathname, ~] = uigetfile(filetype_cells,titleStr,last_browsed_path);
    end
    fullpath = fullfile(pathname, filename);
    if filename ~= 0 % user did not press cancel
        selected_file = fullpath;
        last_browsed_path = pathname;
    else
        selected_file = 0;
    end
end

% Browse for a file and enter it in given edit object
function browse_for_dir(textBox,titleStr)
    global last_browsed_path;
    dirname = uigetdir(last_browsed_path,titleStr); % TODO: better place to start browsing from
    if dirname ~= 0 % user did not press cancel
        set(textBox,'String',dirname);
        last_browsed_path = dirname;
    end
end

% Browse for a file and enter it in given edit object
function copy_dir_from_another_input(destBox,srcBox)
    global last_browsed_path;
    % TODO: if invalid path, not quite sure what this'll do
    srcPath = string_from_box(srcBox);
    [pathPrefix,~,~] = fileparts(srcPath);
    set(destBox,'String',pathPrefix);
    last_browsed_path = pathPrefix; % next browse will be from this location
end

function start_algo(handles)
    try
        % main file reading
        ssa_file = string_from_box(handles.ssa_name);
        ddr_file = string_from_box(handles.ddr_name);
        sb_file = string_from_box(handles.sb_name);
        wa_file = string_from_box(handles.wa_name);
        
        % main file subsetting
        
        % TODO: can we error-check these values at all?
        rows_use_all = bool_from_box(handles.ssa_row_all);
        if rows_use_all
            rows_subset_bounds = [0,0]; % won't be used, doesn't matter
        else
            rows_subset_bounds = [ int_from_box(handles.ssa_row_min), int_from_box(handles.ssa_row_max) ];
        end
        
        cols_use_all = bool_from_box(handles.ssa_col_all);
        if cols_use_all
            cols_subset_bounds = [0,0]; % won't be used, doesn't matter
        else
            cols_subset_bounds = [ int_from_box(handles.ssa_col_min), int_from_box(handles.ssa_col_max) ];
        end
        
        bands_bounds = [ int_from_box(handles.ssa_band_min), int_from_box(handles.ssa_band_max) ];
        
        % general geometry panel
        
        is_L_data = bool_from_box(handles.l_band_button);
        g1_flag = bool_from_box(handles.g1_checkbox);
        
        kernel_size = int_from_box(handles.kernel_size,1,inf); % requires positive value
        pix_spacing = double_from_box(handles.pix_spacing,eps(0),inf); % requires positive value
        pre_end = double_from_box(handles.pre_end,eps(0),inf); % requires positive value
        post_start = double_from_box(handles.post_start,eps(0),inf); % requires positive value
        pix_spacing_output = double_from_box(handles.pix_spacing_output,eps(0),inf); % requires positive value
        init_c_guess = double_from_box(handles.init_ssa_guess,eps(0),inf); % requires positive value
        spec_offset = double_from_box(handles.spec_offset,-inf,inf); % allows any value
        
        proj_start_flag = bool_from_box(handles.projected_start_flag);
        if proj_start_flag
            proj_start_file = string_from_box(handles.proj_scene_name);
        else
            proj_start_file = ''; % won't be used, doesn't matter
        end
        
        % Decomposed penalty panel
        spec_pen_flag = bool_from_box(handles.spec_pen_flag);
        spat_pen_flag = bool_from_box(handles.spat_pen_flag);
        
        beta_spec = double_from_box(handles.beta_spec,eps,inf);
        delta_spec = double_from_box(handles.delta_spec,eps,inf);
        spec_span = int_from_box(handles.spec_span,1,inf);

        beta_spat = double_from_box(handles.beta_spat,eps,inf);
        delta_spat = double_from_box(handles.delta_spat,eps,inf);

        spec_range = double_from_box(handles.spec_range,eps(0),inf);
        pen_internal_iters = int_from_box(handles.pen_internal_iters,1,inf);
        
        % iterations panel
        
        num_iters = int_from_box(handles.num_iters,1,inf);
        switchover_flag = ~bool_from_box(handles.no_switch_radio);
        
        if switchover_flag
            
            % TODO: the switchover-at-iter shouldn't be coupled into penalty
            % switchover in specific like it is in the GUI display
            switch_after_iters = int_from_box(handles.switchover_at_iter);
            if switch_after_iters >= num_iters
                error('Switchover cannot occur after the run is already finished. (Must occur no later than directly after iteration %i)\n',num_iters-1);
            end
            
            % based on which radio button is on, setup that kind of switchover
            if bool_from_box(handles.switch_penalty_radio) % Turn-penalty-on switchover
                switch_new_values = '';
                
                % add whichever kinds of penalty the checkboxes say user wants
                if bool_from_box(handles.switch_penalty_spectral)
                   switch_new_values = strcat(switch_new_values, 'params.spectral_penalty_flag=1;');
                   if(spec_pen_flag) % user also requested spec pen from start: what do they want this for?
                       % TODO: should probably agree on a gentler way to
                       % interpret user. Maybe a default choice of which we
                       % want?
                       error('Both starting spectral penalty and switchover spectral penalty are on. At most one is allowed.');
                   end
                end
                if bool_from_box(handles.switch_penalty_spatial)
                   switch_new_values = strcat(switch_new_values, 'params.spatial_penalty_flag=1;');
                   if(spat_pen_flag) % user also requested spec pen from start: what do they want this for?
                       % TODO: should probably agree on a gentler way to
                       % interpret user. Maybe a default choice of which we
                       % want?
                       error('Both starting spatial penalty and switchover spatial penalty are on. At most one is allowed.');
                   end
                end
                
                % if they didn't ask for anything, probably they meant to
                if(isempty(switch_new_values))
                    error('Switchover for turning on penalty selected, but no types of penalty to turn on were chosen');
                end
                
            else
                error('Unknown type of switchover selected. Most likely this was caused by maintainer error.\n');
            end
            
        else % bool_from_box(handles.no_switch_radio) must be on
            switch_after_iters = 0;
            switch_new_values = '';
        end
        
        % Save data panel
        
        output_prefix = string_from_box(handles.out_file_prefix);
        save_c_flag = bool_from_box(handles.save_c_flag);
        save_mu_flag = bool_from_box(handles.save_mu_flag);
        save_proj_flag = bool_from_box(handles.save_proj_flag);
        save_input_flag = bool_from_box(handles.save_input_flag);
        save_spectrogram_flag = bool_from_box(handles.save_spectrogram_flag);
        if save_c_flag || save_mu_flag || save_proj_flag
            % acquire the non-neg ints from the expression (regardless of
            % delimiters) and make a vector out of them.
            
            raw_save_iters_text = string_from_box(handles.save_iters);
            
            matches = regexp(raw_save_iters_text,'(\d+)','match');
            
            if isempty(matches)
                error('Data saving is on, but no iterations are provided to save on\n');
            end
            
            save_iters = str2double(matches);
            
            % TODO: check that these iters are in range (and whether to
            % error or warn if they're not)
            
        else
            save_iters  = [];
        end
        output_dir = string_from_box(handles.output_dir);
        
        % Geometry Selection panel
        
        off_nadir_flag = bool_from_box(handles.off_nadir_flat_geo_button) || bool_from_box(handles.off_nadir_rolling_geo_button);
        if off_nadir_flag
            radius_file = string_from_box(handles.radius_name);
            topo_file = string_from_box(handles.topo_name);
            spice_file = string_from_box(handles.spice_name);
        else
            radius_file = '';
            topo_file = '';
            spice_file = '';
        end
        
        % Flatness: only not flat if rolling off-nadir geometry selected
        flat = ~bool_from_box(handles.off_nadir_rolling_geo_button);

        % The non-user controlled params
        precompute = false;
        
        % determine whether to use MEX C++ translations & tell user about it
        [use_mex_versions,message] = has_mex();
        if use_mex_versions
            fprintf('C++ MEX translations will be used (they appear to be set up right).\n');
        else
            fprintf('C++ MEX translations will NOT be used. Reason: %s\n', message);
        end
        
        input_start_in = ''; % should be fine if they used absolute paths to input files (TODO)

        disp('Submitting job...');
        
    catch ERR
        % TODO: We should really write errors to the GUI, or pop up a
        % dialog or something. That would probably be easier on user
        fprintf('Error starting run: %s\n',ERR.message);
        return; % no sense attempting to run the main function now
    end
    
    % %%%%%%%%%%%%%%%%%%%%%
    % % STARTING THE RUN! %
    % %%%%%%%%%%%%%%%%%%%%%
    
    % The function call from the black lagoon
        
    % to be safe, we're giving explicit values for ALL parameters the
    % main code recognizes (no tricky-to-find defaults for us!)
    
    test5v5_compact(...
        'spatial_penalty_flag', spat_pen_flag,...
        'spectral_penalty_flag', spec_pen_flag,...
        'beta_spat', beta_spat,...
        'delta_spat', delta_spat,...
        'beta_spec', beta_spec,...
        'delta_spec', delta_spec,...
        'spectral_span', spec_span,...
        'max_iter', pen_internal_iters,...
        'is_L_data', is_L_data,...
        'size_kernel', kernel_size,...
        'pix_spacing', pix_spacing,...
        'pre_end', pre_end,...
        'post_start', post_start,...
        'spec_range',spec_range,...
        'pix_spacing_output', pix_spacing_output,...
        'pix_value', init_c_guess,...
        'num_iter', num_iters,...
        'g1_flag', g1_flag,...
        'geoflags.offnadir', off_nadir_flag,...
        'geoflags.precompute', precompute,...
        'geoflags.flat', flat,...
        'use_mex_versions', use_mex_versions,...
        'switchover_flag', switchover_flag,...
        'switchover_at_iter', switch_after_iters,...
        'switchover_code', switch_new_values,...
        'spectral_offset', spec_offset,...
        'save_c_on_iters_flag', save_c_flag,...
        'save_mu_on_iters_flag', save_mu_flag,...
        'save_proj_on_iters_flag', save_proj_flag,...
        'save_input_flag', save_input_flag,...
        'save_spectrogram_flag', save_spectrogram_flag,...
        'saving_iters', save_iters,...
        'crism_iof_filename', ssa_file,...
        'ddr_iof_filename', ddr_file,...
        'sb_filename', sb_file,...
        'wa_filename', wa_file,...
        'filenames.radius_filename', radius_file,...
        'filenames.topo_filename', topo_file,...
        'filenames.spice_filename', spice_file,...
        'start_from_scene_flag', proj_start_flag,...
        'start_from_scene_filename', proj_start_file,...
        'input_path_start_in', input_start_in,...
        'output_path_start_in', output_dir,...
        'rows_use_all', rows_use_all,...
        'rows_subset_bounds', rows_subset_bounds,...
        'cols_use_all', cols_use_all,...
        'cols_subset_bounds', cols_subset_bounds,...
        'bands_bounds', bands_bounds,...
        'observation_base', output_prefix...
        );

end

function auto_fill_sb_wa(handles)
% attempt to auto-fill sb & wa fields, if those files in same dir as ssa file
    ssaFullPath = string_from_box(handles.ssa_name);
    [ssaDir,~,~] = fileparts(ssaFullPath);
    
    % what we're looking for
    is_L_data = bool_from_box(handles.l_band_button);
    if is_L_data
        sb_filename = 'CDR490947778566_SB0000000L_3.IMG';
        wa_filename = 'CDR490947778566_WA0000000L_3.IMG';
    else
        sb_filename = 'CDR450924300802_SB0000000S_2.IMG';
        wa_filename = 'CDR450924300802_WA0000000S_2.IMG';
    end
    
    found_sb = false;
    found_wa = false;
    
    % go over possibilities
    files_in_dir = dir(ssaDir);
    for i=1:length(files_in_dir)
        
        next_file = files_in_dir(i).name;
        
        if strcmpi(next_file,sb_filename)
            found_sb = true;
            if found_wa
                break; % both found: we're done here
            end
        end
        
        if strcmpi(next_file,wa_filename)
            found_wa = true;
            if found_sb
                break; % both found: we're done here
            end
        end
        
    end
    
    if found_sb && found_wa % They're in there
        set(handles.sb_name,'String',fullfile(ssaDir,sb_filename));
        set(handles.wa_name,'String',fullfile(ssaDir,wa_filename));
    else
       % Could not find files. Design decision: just print msg to console
       fprintf('Warning: The corresponding SB & WA files couldn''t be found in the directory of the SSA file.\n');
    end
    
end

function save_options_file(handles, opts_filename)
    opts_str = '';
    handle_names = fieldnames(handles);
    for i = 1:numel(handle_names)
        current_handle = handles.(handle_names{i});
        try
            current_style = get(current_handle,'Style');
        catch
            current_style = ''; % something where style is not a valid property: we don't care about it then
        end
        
        switch current_style % Type of box
            case 'edit' % textbox to type in
                opts_str = [opts_str sprintf('%s:%s\n',handle_names{i},string_from_box(current_handle))];
            case 'checkbox'
                opts_str = [opts_str sprintf('%s:%i\n',handle_names{i},bool_from_box(current_handle))];
            case 'radiobutton'
                opts_str = [opts_str sprintf('%s:%i\n',handle_names{i},bool_from_box(current_handle))];
            otherwise
                % other types are static text, etc and have nothing
                % interesting to store
        end
    end
    
    fid = fopen(opts_filename,'w+');
    fprintf(fid,'%s',opts_str);
    fclose(fid);
    
end

function restore_options_file(handles, opts_filename)
% read options from the given file. Expect exactly 1 option saved per line
    fid = fopen(opts_filename);
    fprintf('Loading GUI options from %s...\n',opts_filename);
    next_line = fgetl(fid);
    while ~(isnumeric(next_line) && next_line == -1) % file has lines left
        result = regexp(next_line,'(\w+):(.*)$','tokens');
        if ~isempty(result) % there is a match in this line
            handle_name = result{1}{1};
            value = result{1}{2};
            if isfield(handles,handle_name) % the GUI has a field by this name: set it.
            current_handle = handles.(handle_name);
            switch current_handle.Style
                case 'edit'
                    set(current_handle,'String',value);
                    gui_react(current_handle,handles); % to maintain consistant GUI state (in terms of opaqueness of boxes, etc)
                case 'checkbox'
                    val = str2double(value);
                    if isnan(val)
                        val = 0;
                    end
                    set(current_handle,'Value',val);
                    gui_react(current_handle,handles); % to maintain consistant GUI state (in terms of opaqueness of boxes, etc)
                case 'radiobutton'
                    val = str2double(value);
                    if isnan(val)
                        val = 0;
                    end
                    set(current_handle,'Value',val);
                    gui_react(current_handle,handles); % to maintain consistant GUI state (in terms of opaqueness of boxes, etc)
                otherwise
                    fclose(fid);
                    error('The field %s contained in the options file cannot take user input\n');
            end
            else
                % the GUI has no field: warn
                fprintf('Warning: ignoring field ''%s'', which is not supported by the GUI. This option-file\nmight be for a previous version of this program.\n',handle_name);
            end
        end
        next_line = fgetl(fid);
    end
    fprintf('Loading GUI options from %s... done.\n',opts_filename);
    fclose(fid);
end

function s = string_from_box(h)
    s = get(h,'String');
end

function d = double_from_box(h,rangeMin,rangeMax)
    % [rangeMin, rangeMax] forms closed interval of allowed values. If not
    % given, the min defaults to 0 and the max to Inf (all positive values)
    if nargin < 3
        if nargin < 2
            rangeMin = 0;
        end
        rangeMax = Inf;
    end
    
    d = str2double(string_from_box(h));
    if isnan(d)
        % String was not a number
        error('The form box %s must contain a number\n',get(h,'Tag'));
    end
    
    if d>rangeMax || d<rangeMin
        error('The form box %s requires a value in the range %f to %f\n',get(h,'Tag'),rangeMin,rangeMax);
    end
end

function i = int_from_box(h,rangeMin,rangeMax)
    % [rangeMin, rangeMax] forms closed interval of allowed values. If not
    % given, the min defaults to 0 and the max to Inf (all positive values)
    if nargin < 3
        if nargin < 2
            rangeMin = 0;
        end
        rangeMax = Inf;
    end
    
    % Get the value as a double
    i = double_from_box(h,rangeMin,rangeMax);
    
    % Make sure value is integral
    if i ~= floor(i)
        error('The form box %s requires an integer value.\n',get(h,'Tag'));
    end
end

function b = bool_from_box(h)
    b = get(h,'Value');
end

% To get the vectors or 'all' that are currently used to index main file
% dimensions
function vect = get_formatted_int_range(h_all_button,h_min,h_max)
    if bool_from_box(h_all_button)
        vect = 'all';
    else
        vect = double_from_box(h_min):double_from_box(h_max);
    end
end
