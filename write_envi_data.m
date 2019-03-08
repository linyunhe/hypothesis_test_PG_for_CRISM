function write_envi_data(data,base_filename,interleave,varargin)
%WRITE_ENVI_DATA writes file with data and a corresponding ENVI header
% Takes a matrix of data, a filename & interleave for the data, and maybe
% additional arguments to be passed to save_envi_header (giving information
% like description text or wavelengths to be saved in the header). See that
% function's notes for a listing of understood options.
    
    % Get filename for header
    [path,mainname,ext] = fileparts(base_filename);
    if strcmp(ext,'')
        % if no extension, build generic ones for data & header
        data_filename = strcat(base_filename, '.', interleave);
        header_filename = strcat(base_filename, '.hdr');
    else
        % if there is an extension, replace it for header file
        data_filename = base_filename;
        header_filename = fullfile(path,sprintf('%s.hdr',mainname));
    end
    
    args = parseParameters(varargin);
    
    % write data
    if isfield(args, 'type')
        % write in requested precision
        multibandwrite(data, data_filename, interleave, 'precision', args.type);
    else
        % write in precision the data came in
        multibandwrite(data, data_filename, interleave);
    end

    % write header
    save_envi_header(data, header_filename, interleave, varargin);

end
