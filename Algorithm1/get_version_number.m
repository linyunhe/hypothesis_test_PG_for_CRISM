function [vernum,veralgo] = get_version_number()
% retrieves the number of this version of CRISM (e.g. 4.4.0) and the algorithm type (e.g. MLM), stored in a local file
    fid = fopen('VERSION');
    if fid ~= -1 % opened alright
        vernum = fgetl(fid);
        veralgo = fgetl(fid);
        fclose(fid);
    else
        vernum = '<VER>';
        veralgo = '<ALGO>';
    end
end