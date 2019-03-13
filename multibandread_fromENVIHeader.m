function out_matrix = multibandread_fromENVIHeader( data_filename )
%multibandread_fromENVIHeader Reads a multiband datafile using info from its ENVI header file
%   Detailed explanation goes here

% Get the header data. Header file name is this one's name w/ extension
% replaced
header_filename = get_header_filename(data_filename);
header_data = read_envi_header(header_filename);

data_dims = [header_data.lines,header_data.samples,header_data.bands];

% Test to make sure that the file really contains the amount of data
% its header claims
%validate_datafile_size(data_filename, data_dims, header_data.dataType);

% NOTE: here assuming offset of 0 and little-endian, because this data has
% always been like that. This function and read_ENVI_header will require
% adjustment if that ever changes

out_matrix = multibandread(data_filename,data_dims,...
    header_data.dataType,0,header_data.interleave,'ieee-le');

end

function header_filename = get_header_filename(data_filename)
    [pathstr,mainname,ext] = fileparts(data_filename);
    
    if isempty(pathstr)
        pathstr = '.';
    end
    
    hypothetical_header_filename = sprintf('%s.hdr',mainname);
    
    % our guess at the filename might have incorrect case (and this matters
    % in Linux). See if any of the data file's siblings are of the correct
    % form.
    
    directory_contents = getContentsList(pathstr);
    
    files_match_vect = strcmpi(hypothetical_header_filename,directory_contents); % case-insensitive
    
    % error checking
    if sum(files_match_vect)>1
        error('There are multiple candidates for the header file %s in path %s',hypothetical_header_filename,pathstr);
    elseif sum(files_match_vect)==0
        error('There are no candidates for the header file %s in path %s',hypothetical_header_filename,pathstr);
    end
    
    header_filename = fullfile(pathstr,directory_contents{files_match_vect});
    
end

function filenames = getContentsList(pathstr)
    fileStruct = dir(pathstr);
    filenames = cell(length(fileStruct),1);
    for i=1:length(fileStruct)
        filenames{i} = fileStruct(i).name;
    end
end