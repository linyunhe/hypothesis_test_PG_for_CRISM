function mro_position = spice_info(SPICE_filename, first_row, last_row)
%Read in SPICE data, calculated from WebGeoCalc, which gives the
%spacecraft's position in time increments chosen to match the time each
%line of the image was taken

    fid = fopen(SPICE_filename);
    %Skip the header lines
    for headerline = 1:21
        fgetl(fid);
    end
    %Read data as strings into cell array
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s';
    celldata = textscan(fid,formatSpec);
    fclose(fid);

    lon_mro = str2num(char(celldata{4}))*pi/180;
    phi_mro = lon_mro;
    lat_mro = str2num(char(celldata{5}))*pi/180;
    theta_mro = pi/2 - lat_mro;
    r_mro = str2num(char(celldata{6}))*1000; %radius of MRO from center of Mars
    
    mro_position = zeros((last_row-first_row+1), 3);
    %read in only the desired rows
    for row = first_row:last_row
        mro_position(row-first_row+1,1) = r_mro(row)*sin(theta_mro(row))*cos(phi_mro(row));
        mro_position(row-first_row+1,2) = r_mro(row)*sin(theta_mro(row))*sin(phi_mro(row));
        mro_position(row-first_row+1,3) = r_mro(row)*cos(theta_mro(row));
    end
end