function validate_datafile_size( filepath, expecteddimensions, precision )
%VALIDATE_DATAFILE_SIZE Stops program if datafile has unexpected size
%   Terminates program if the given file does not have the expected
%   length of elements of the given type.
%
%   Author: Daniel Politte (9 June 2015)
    fid = fopen(filepath,'r');
    [~,count] = fread(fid,inf,precision);
    fclose(fid);
    
    if count ~= prod(expecteddimensions)
        error('Data file has unexpected size. Check the dimension parameters.');
    end
end

