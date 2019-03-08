function save_envi_header_mapproject(data,filename,interleave,desc_text,wavelengths,mapProjectionInfo)
    % wavelengths is a vector with length of the # of bands in the data
    % (dimension 3), ordered starting with band #1.
    fid = fopen(filename,'w');
    
    [lines,samples,bands] = size(data);
    
    % error conditions
    if bands~=length(wavelengths)
        disp('error in save_envi_header: mismatch between size of data and given wavelengths');
        return;
    end
    
    if ~any(strcmp(interleave,{'bsq','bil','bip'}))
        fprintf('error in save_envi_header: invalid interleaving method %s\r\n',interleave);
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
    
    % Map projection info strings
    mapInfoString = sprintf('map info = {Mars Equirectangular Default, %f, %f, %f, %f, %f, %f, MARS SphereD, units=Meters, rotation=%f}', ...
        mapProjectionInfo.tieRow, mapProjectionInfo.tieCol, ...
        mapProjectionInfo.tieLonM, mapProjectionInfo.tieLatM, ...
        mapProjectionInfo.pixelSizeX, mapProjectionInfo.pixelSizeY, ...
        mapProjectionInfo.rotation);
    projInfoString = 'projection info = {17, 3396190.0, 0.000000, 0.000000, 0.0, 0.0, MARS SphereD, Mars Equirectangular Default, units=Meters}';
    coSysString = 'coordinate system string = {PROJCS["Mars Equirectangular Default",GEOGCS["GCS_Unknown",DATUM["D_Unknown",SPHEROID["S_Unknown",3396190.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Equidistant_Cylindrical"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",0.0],UNIT["Meter",1.0]]}';
    
    
    wavelengths_str = wavelengths_str(1:length(wavelengths_str)-2); % kill trailing comma-space
    
    % print out to file
    
    fprintf(fid,'ENVI\r\ndescription = {\r\n%s\r\n}\r\n', desc_text);
    fprintf(fid,'samples = %i\r\nlines   = %i\r\nbands   = %i\r\n',samples,lines,bands);
    fprintf(fid,'header offset = 0\r\n');
    fprintf(fid,'file type = ENVI Standard\r\n');
    fprintf(fid,'data type = 5\r\n');
    fprintf(fid,'interleave = %s\r\n',interleave);
    fprintf(fid,'sensor type = Unknown\r\n');
    fprintf(fid,'byte order = 0\r\n');
    fprintf(fid,'wavelength units = Micrometers\r\n');
    fprintf(fid,'%s\r\n',mapInfoString);
    fprintf(fid,'%s\r\n',projInfoString);
    fprintf(fid,'%s\r\n',coSysString);
    fprintf(fid,'wavelength = {');

    fprintf(fid,'%s',wavelengths_str);
    
    fprintf(fid,'}');
    fprintf(fid,'\r\n');
    
    fclose(fid);
end