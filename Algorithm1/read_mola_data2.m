function scene_areoid = read_mola_data2(filenames,c_num_rows,c_num_cols,ddr_lat,ddr_lon)

%Make a topography map
radius_label_filename = strcat([filenames.radius_filename(1:end-3),'lbl']);
topo_label_filename = strcat([filenames.topo_filename(1:end-3),'lbl']);

scene_radius = read_file(filenames.radius_filename,radius_label_filename,c_num_rows,c_num_cols,ddr_lat,ddr_lon);
scene_topo = read_file(filenames.topo_filename,topo_label_filename,c_num_rows,c_num_cols,ddr_lat,ddr_lon);

%substract them to get the baseline
scene_areoid = scene_radius-scene_topo;

end

function projection = read_file(filename,label_filename,c_num_rows,c_num_cols,ddr_lat,ddr_lon)
%read in data
fid = fopen(filename);
data = fread(fid, 'int16', 'ieee-be');
fclose(fid);
%now interpret it
label_data = fileread(label_filename);
RES = find_value(label_data, 'MAP_RESOLUTION');
SAMPLE_PROJECTION_OFFSET = find_value(label_data, 'SAMPLE_PROJECTION_OFFSET');
CENTER_LONGITUDE = find_value(label_data, 'CENTER_LONGITUDE');
LINE_PROJECTION_OFFSET = find_value(label_data, 'LINE_PROJECTION_OFFSET');
CENTER_LATITUDE = find_value(label_data, 'CENTER_LATITUDE');
LINES = find_value(label_data, 'LINES');
LINE_SAMPLES = find_value(label_data, 'LINE_SAMPLES');
OFFSET = find_value(label_data, 'OFFSET');

data = reshape(data, [LINE_SAMPLES, LINES])';
%create array the size of c to hold MOLA data
projection = zeros(c_num_rows, c_num_cols);
for c_row = 1:c_num_rows
    LAT = ddr_lat(c_row)*180/pi;
    %find the corresponding row in the MOLA arrays from the equation (dsmap.cat)
    LINE = LINE_PROJECTION_OFFSET - RES * (LAT - CENTER_LATITUDE);
    for c_col = 1:c_num_cols
        LON = ddr_lon(c_col)*180/pi;
        if LON < 0
            LON = LON+360;
        end
        %find the corresponding columin in MOLA arrays
        SAMPLE = SAMPLE_PROJECTION_OFFSET + RES * (LON - CENTER_LONGITUDE);
        projection(c_row, c_col) = data(round(LINE), round(SAMPLE)) + OFFSET;
    end
end

end

function value = find_value(label_str,varname)
    equation = strcat([varname,'\s*=\s*-?\d+.?\d*']);
    [startEquation,endEquation] = regexp(label_str,equation);
    subString = label_str(startEquation:endEquation);
    rightSide = '-?\d+.?\d*';
    [startNum,endNum] = regexp(subString,rightSide);
    value = str2num(char(subString(startNum:endNum)));
end