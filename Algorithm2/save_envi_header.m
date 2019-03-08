function save_envi_header(data,filename,interleave,desc_text,wavelengths)
    % wavelengths is a vector with length of the # of bands in the data
    % (dimension 3), ordered starting with band #1.
    fid = fopen(filename,'w');
    
    [lines,samples,bands] = size(data);
    if wavelengths ~= 0
    
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
    
    wavelengths_str = wavelengths_str(1:length(wavelengths_str)-2); % kill trailing comma-space
    else
     
    wavelengths_str = 'Not bands';
    end
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
    fprintf(fid,'wavelength = {');

    fprintf(fid,'%s',wavelengths_str);
    
    fprintf(fid,'}');
    fprintf(fid,'\r\n');
    
    fclose(fid);
end