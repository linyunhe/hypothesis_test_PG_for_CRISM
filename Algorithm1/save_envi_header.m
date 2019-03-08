function save_envi_header(data,filename,interleave,varargin)
    % wavelengths is a vector with length of the # of bands in the data
    % (dimension 3), ordered starting with band #1.
    % TODO: document varargs stuff more clearly
    
    % Possible parameters in varargin:
    % - 'desc': text to use, if not there then generic msg
    % - 'type': the type to translate into ENVI format. If absent, than use
    % the type of data param
    % - 'wavelengths': a vector to use. If absent, omit section altogether.
    % - 'default bands': an RGB vector to use (length 3). If absent, omit.
    % - 'mapProjectionInfo': an object describing map projection
    % information relavent to the saved data. If absent, omit.
    % - 'badbandslist': a vector to put as the bad bands list (bbl). If
    % absent, omit section.
    % - 'band_names': a list of names for the bands given as a cell array.
    % Maps into the header record type labeled "band names".
    
    % handle extra args being passed as a single cell array rather than multiple args
    if iscell(varargin) 
        % TODO: what if there's a cell array _and_ some proper args? This will disregard all but the first item then.
        varargin = varargin{1}; 
    end
    
    % default vals and things
    desc = 'ENVI header saved using CRISM MLM or related code';
    type_str = class(data); % whatever type the data is
    wavelengths = [];
    default_bands = [];
    mapProjectionInfo = [];
    badbands = [];
    bandNames = {};
    
    % fill over defaults using given values
    nargs = length(varargin);

    for i=1:2:nargs-1
        % We expect that a parameter name is at i, and the value is at i+1

        if ~ischar(varargin{i}) % if paramter name isn't a string, abort
            error('Expected argument #%i to be a string (command field name): instead it''s a %s',...
                i,class(varargin{i}));
        end
        
        switch varargin{i}
            case 'desc'
                desc = varargin{i+1};
            case 'type'
                type_str = varargin{i+1};
            case 'wavelengths'
                wavelengths = varargin{i+1};
            case 'default_bands'
                default_bands = varargin{i+1};
                if length(default_bands) ~= 3
                    if isempty(default_bands)
                        % TODO: this handling is suboptimal
                        %disp('Warning: default bands not printed in output file');
                    else
                        error('Vector of default bands has length %i (should be 3, with entries being RGB respectively)',length(default_bands));
                    end
                end
            case 'mapProjectionInfo'
                mapProjectionInfo = varargin{i+1};
                % TODO: validate that struct has relavent properties
            case 'badbandslist'
                badbands = varargin{i+1};
            case 'band_names'
                bandNames = varargin{i+1};
            otherwise
                error('Unrecognized argument ''%s''',varargin{i});
        end

    end
    
    % translate given type string to ENVI's id number
    switch type_str
        case 'single'
            type = 4; % not 'single', because these types meant for multibandread
        case 'double'
            type = 5;
        otherwise
            error('There is no know ENVI datatype code for the type ''%s''. If this is a valid type, please update this function.', type_str);
    end
    
    % error conditions
    if ~any(strcmp(interleave,{'bsq','bil','bip'}))
        error('error in save_envi_header: invalid interleaving method ''%s''',interleave);
    end

    [lines,samples,bands] = size(data);
    
    if ~isempty(wavelengths) % only if call specified some wavelengths to use
        % error condition
        if bands~=length(wavelengths)
            disp('error in save_envi_header: mismatch between size of data and given wavelengths');
            return;
        end

        wavelengths = wavelengths / 1000; % convert from nm to um

        wavelengths_per_line = 7;

        % build a string of the wavelength values, comma separated, with a
        % reasonable number per line
        wavelengths_str = '';
        for line_start_index = 1:wavelengths_per_line:bands
            new_wavelengths = sprintf('%f, ',wavelengths(line_start_index:min(line_start_index + wavelengths_per_line - 1, bands)));
            wavelengths_str = sprintf('%s\r\n%s',wavelengths_str,new_wavelengths);
        end

        wavelengths_str = wavelengths_str(1:length(wavelengths_str)-2); % kill trailing comma-space
    end
    
    if ~isempty(mapProjectionInfo) % only if call specified map projection to use
        % Map projection info strings
        mapInfoString = sprintf('map info = {Mars Equirectangular Default, %f, %f, %f, %f, %f, %f, MARS SphereD, units=Meters, rotation=%f}', ...
            mapProjectionInfo.tieRow, mapProjectionInfo.tieCol, ...
            mapProjectionInfo.tieLonM, mapProjectionInfo.tieLatM, ...
            mapProjectionInfo.pixelSizeX, mapProjectionInfo.pixelSizeY, ...
            mapProjectionInfo.rotation);
        projInfoString = 'projection info = {17, 3396190.0, 0.000000, 0.000000, 0.0, 0.0, MARS SphereD, Mars Equirectangular Default, units=Meters}';
        coSysString = 'coordinate system string = {PROJCS["Mars Equirectangular Default",GEOGCS["GCS_Unknown",DATUM["D_Unknown",SPHEROID["S_Unknown",3396190.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Equidistant_Cylindrical"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",0.0],UNIT["Meter",1.0]]}';
    end
    
    % print out to file
    
    fid = fopen(filename,'w');
    
    fprintf(fid,'ENVI\r\ndescription = {\r\n%s\r\n}\r\n', desc);
    fprintf(fid,'samples = %i\r\nlines   = %i\r\nbands   = %i\r\n',samples,lines,bands);
    fprintf(fid,'header offset = 0\r\n');
    fprintf(fid,'file type = ENVI Standard\r\n');
    fprintf(fid,'data type = %i\r\n',type);
    fprintf(fid,'interleave = %s\r\n',interleave);
    fprintf(fid,'sensor type = Unknown\r\n');
    fprintf(fid,'byte order = 0\r\n');
    if ~isempty(wavelengths)
        fprintf(fid,'wavelength units = Micrometers\r\n');
        fprintf(fid,'wavelength = {');

        fprintf(fid,'%s',wavelengths_str);

        fprintf(fid,'}\r\n');
    end
    
    if ~isempty(default_bands)
        fprintf(fid,'default bands={%i,%i,%i}\r\n',default_bands(1),default_bands(2),default_bands(3));
    end
    
    if ~isempty(mapProjectionInfo)
        fprintf(fid,'%s\r\n',mapInfoString);
        fprintf(fid,'%s\r\n',projInfoString);
        fprintf(fid,'%s\r\n',coSysString);
    end
    
    if ~isempty(badbands)
        fprintf(fid,'bbl={\n'); % start
        bandsstr = sprintf('%i,',badbands); % all the bands
        bandsstr = bandsstr(1:end-1); % cut off trailing comma
        fprintf(fid,bandsstr); % printing the bands
        fprintf(fid,'\n}\n'); % end
    end
    
    if ~isempty(bandNames)
        bandNamesStr = 'band names={\n';
        for bn = bandNames
            bandNamesStr = [bandNamesStr, sprintf('%s,',bn{1})];
        end
        bandNamesStr = bandNamesStr(1:end-1); % kill trailing comma
        bandNamesStr = [bandNamesStr, sprintf('}\n')];
        fprintf(fid,bandNamesStr);
    end
    
    fprintf(fid,'\r\n');
    fclose(fid);
end
