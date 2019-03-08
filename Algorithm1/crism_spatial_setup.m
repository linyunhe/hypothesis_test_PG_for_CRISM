function [c,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_positions,pix_normals,output_mask] = crism_spatial_setup(geoflags,filenames,ddr_iof,pix_spacing,pix_value,num_rows,num_cols,num_bands,kernel_size)

%Projects sensor space onto the surface and interpolates to make an initial
%guess at a scene (all constant values, equal to pix_value). Also generates
%useful mappings between sensor space pixels and scene pixels. If geometry
%is included, also calls another function to generate (and maybe
%interpolate) a topography map.

r = 3396190.0; %radius of Mars in meters (exact value used in SuperGLT reconstruction)
span = (kernel_size-1)/2; % for filling C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate dimensions of c image%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%min/max lat in mu
min_lat = min(min(ddr_iof(:,:,4)*pi/180.0));
max_lat = max(max(ddr_iof(:,:,4)*pi/180.0));

%min/max lon in mu
min_lon = min(min(ddr_iof(:,:,5)*pi/180.0));
max_lon = max(max(ddr_iof(:,:,5)*pi/180.0));

%vertical distance in meters between min and max lat
y_dist = r*(max_lat-min_lat);
%horizotnal distance in meters between min and max lon
x_dist = r*sin(pi/2-abs(max_lat+min_lat)/2)*(max_lon-min_lon);

%find the number of rows and columns in c based on desired pixel size
c_num_rows = 1+round(y_dist / pix_spacing);
c_num_cols = 1+round(x_dist / pix_spacing);

c = zeros(c_num_rows, c_num_cols, num_bands); %fill regions outside bow tie with default value 0.0 
% mask = c;
%find the difference in latitude between neighboring rows and the difference in longitude between neighboring columns 
lat_spacing = abs(max_lat-min_lat)/(c_num_rows-1); %c_lat_spacing
lon_spacing = abs(max_lon-min_lon)/(c_num_cols-1); %c_lat_spacing

%initialize arrays
lat_number = zeros(num_rows,num_cols); %c_row_index
lon_number = zeros(num_rows,num_cols); %c_col_index

%select top left pixel for indexing
center_row = 1; %J: why is it called center, and why would we ever need to adjust this? Origin in coordinate system?
center_col = 1;

%intialize arrays. Every c row will have the same latitude and every c
%column will have the same longitude.

%determine latitude for every c row given that c rows are evenly spaced
ddr_lat = min_lat+lat_spacing*((1:c_num_rows)-center_row);

%determine longitude for every c col given that c columns are evenely spaced
ddr_lon = max_lon-lon_spacing*((1:c_num_cols)-center_col);

%set pixels inside of the bow tie to 0.1. This can be any value greater
%than 0 and less than 1.
for mu_row=1:num_rows
    for mu_col=1:num_cols
        
        %determine row in c that mu pixel corresopnds to
        lat_ratio = abs(ddr_iof(mu_row,mu_col,4)*pi/180.0-ddr_lat(1))/lat_spacing; %should this 1 really be center_row? Also rename lat_ratio as another kind of index in c
        lat_number(mu_row,mu_col) = 1+ round(lat_ratio); %where does the extra 1 come from?
        
        %determine col in c that mu pixel corresponds to
        lon_ratio = abs(ddr_iof(mu_row,mu_col,5)*pi/180.0-ddr_lon(1))/lon_spacing;
        lon_number(mu_row,mu_col) = 1+ round(lon_ratio);
        
        % fill in pixels in C around this pixel's location equivalent to
        % the size of the spatial kernel we'll use
        row_min = max(1,lat_number(mu_row,mu_col)-span);
        row_max = min(lat_number(mu_row,mu_col)+span,c_num_rows);
        col_min = max(1,lon_number(mu_row,mu_col)-span);
        col_max = min(lon_number(mu_row,mu_col)+span,c_num_cols);
        c(row_min:row_max,col_min:col_max,:) = pix_value;

    end
end

% Setup a mask for projected-scene output. The mask should include no
% edge-pixels whose values have a lot of artifacts.
se = strel('square',3);
output_mask = c(:,:,1)>0;

for i=1:span
    output_mask = imerode(output_mask,se);
end

% Make sure all edge pixels removed (while sacrificing a small number of
% good ones)
output_mask(1:span,:) = 0;
output_mask(end-span+1:end,:) = 0;
output_mask(:,1:span) = 0;
output_mask(:,end-span+1:end) = 0;


c_temp=c(:,:,1); %Specify why we need this

mapping = zeros(num_rows,2);
mapping(:,1) = ddr_iof(:,1,4)*pi/180.0;

for i=1:num_rows
    mapping(i,2)=i;
end

order = sortrows(mapping, 1);

%If using new geometry, find the radius of the areoid (~sea level) in this
%area on Mars, using the MOLA MEGDRs
if geoflags.offnadir
    scene_areoid = read_mola_data2(filenames,c_num_rows,c_num_cols,ddr_lat,ddr_lon); %a big array
    scene_areoid = mean(mean(scene_areoid)); %Since the MEGDR is low-resolution, we take the mean to avoid big discontinuities in the scene
end

for n=1:num_cols-1,
    for i=1:num_rows-1,
       lat_new_diff = abs(order(i,1)-order(i+1,1))/4.000;
       lon_new_diff = abs(ddr_iof(order(i,2),n,5)*pi/180.0-ddr_iof(order(i+1,2),n,5)*pi/180.0)/4.000;
       for t=1:4
           fill_c_row = lat_number(i,n)+round(lat_new_diff*t/lat_spacing);
           fill_c_col = lon_number(i,n)+round(lon_new_diff*t/lon_spacing);
           c_temp(fill_c_row,fill_c_col,:) = pix_value;
       end
    end
end



%interpolate
for iter=1:8
    for c_row=1:(c_num_rows)
        for c_col=1:(c_num_cols)

            if c_temp(c_row,c_col)~=pix_value

                %boundary conditions
                %%%%%%
                %booleans
                upper = c_row-4<1;               %near the top
                lower = c_row+4>c_num_rows;     %near the bottom
                right = c_col+4>c_num_cols;     %near the right edge
                left = c_col-4<1;               %near the left edge
                
                if ~right
                    right_neighbors = sum(c_temp(c_row,c_col+1:c_col+4))>pix_value; %There is >=1 neighbor in the nearest four right pixels
                end
                
                if ~lower
                    lower_neighbors = sum(c_temp(c_row+1:c_row+4,c_col))>pix_value; %neighbor in lower pixels
                end
                
                if ~left
                    left_neighbors = sum(c_temp(c_row,c_col-4:c_col-1))>pix_value;  %neighbor in left pixels
                end
                
                if ~upper
                    upper_neighbors = sum(c_temp(c_row-4:c_row-1,c_col))>pix_value; %neighbor in upper pixels
                end
                %%%%%%
                
                if upper && left                                            %(c_row-4)<1 && (c_col-4)<1          
                    if right_neighbors && lower_neighbors                   %sum(c_temp(c_row,c_col+1:c_col+4))>pix_value && sum(c_temp(c_row+1:c_row+4,c_col))>pix_value
                        c(c_row,c_col,:) = pix_value;                        
                    end
                    
                elseif upper && right                                       %(c_row-4)<1 && (c_col+4)>c_num_cols                   
                    if left_neighbors && lower_neighbors                    %sum(c_temp(c_row,c_col-4:c_col-1)) > pix_value && sum(c_temp(c_row+1:c_row+4,c_col)) > pix_value
                        c(c_row,c_col,:) = pix_value;                        
                    end
                
                elseif lower && right                                       %c_row+4>c_num_rows && c_col+4>c_num_cols                   
                    if upper_neighbors && left_neighbors                    %sum(c_temp(c_row-4:c_row-1,c_col))>pix_value && sum(c_temp(c_row,c_col-4:c_col-1))>pix_value
                        c(c_row,c_col,:) = pix_value;                        
                    end
                    
                elseif lower && left                                        %(c_row+4)>c_num_rows && (c_col-4)<1   
                    if upper_neighbors && right_neighbors                   %sum(c_temp(c_row-4:c_row-1,c_col))>pix_value && sum(c_temp(c_row,c_col+1:c_col+4))>pix_value
                        c(c_row,c_col,:) = pix_value;                    
                    end
                
                elseif lower                                                %((c_row+2)>c_num_rows)||((c_row+1)>c_num_rows)||((c_row+3)>c_num_rows)||((c_row+4)>c_num_rows)
                    if upper_neighbors && left_neighbors                    %sum(c_temp(c_row-4:c_row-1,c_col))>pix_value && sum(c_temp(c_row,c_col-4:c_col-1))>pix_value
                        c(c_row,c_col,:) = pix_value;
                    elseif upper_neighbors && right_neighbors               %sum(c_temp(c_row-4:c_row-1,c_col))>pix_value && sum(c_temp(c_row,c_col+1:c_col+4))>pix_value
                        c(c_row,c_col,:) = pix_value;                    
                    end

                elseif right                                                %(c_col+4)>c_num_cols%((c_col+2)>c_num_cols)||((c_col+1)>c_num_cols)||((c_col+3)>c_num_cols)||((c_col+4)>c_num_cols)
                    if upper_neighbors && left_neighbors                    %sum(c_temp(c_row-4:c_row-1,c_col))>pix_value && sum(c_temp(c_row,c_col-4:c_col-1))>pix_value
                        c(c_row,c_col,:) = pix_value;     
                        %THIS NEXT LINE ACTUALLY DOESN'T MAKE SENSE: SHOULD
                        %CHECK LOWER AND LEFT neighbors, NOT LOWER AND UPPER
                    elseif lower_neighbors && upper_neighbors               %sum(c_temp(c_row+1:c_row+4,c_col))>pix_value && sum(c_temp(c_row-4:c_row-1,c_col))>pix_value
                        c(c_row,c_col,:) = pix_value;                        
                    end

                elseif upper                                                %(c_row-4)<1
                    if left_neighbors && lower_neighbors                    %sum(c_temp(c_row,c_col-4:c_col-1))>pix_value && sum(c_temp(c_row+1:c_row+4,c_col))>pix_value
                        c(c_row,c_col,:) = pix_value;
                    elseif right_neighbors && lower_neighbors               %sum(c_temp(c_row,c_col+1:c_col+4))>pix_value && sum(c_temp(c_row+1:c_row+4,c_col))>pix_value
                        c(c_row,c_col,:) = pix_value;                        
                    end

                elseif left                                                 %(c_col-4)<1
                    if upper_neighbors && right_neighbors                   %sum(c_temp(c_row-4:c_row-1,c_col))>pix_value && sum(c_temp(c_row,c_col+1:c_col+4))>pix_value
                        c(c_row,c_col,:) = pix_value;                        
                    elseif right_neighbors && lower_neighbors               %sum(c_temp(c_row,c_col+1:c_col+4))>pix_value && sum(c_temp(c_row+1:c_row+4,c_temp))>pix_value
                        c(c_row,c_col,:) = pix_value;                        
                    end

                    
                %actual interpolation
                    
                elseif upper_neighbors && left_neighbors                    %sum(c_temp(c_row-4:c_row-1,c_col))>pix_value && sum(c_temp(c_row,c_col-4:c_col-1))>pix_value
                    c(c_row,c_col,:) = pix_value;

                elseif upper_neighbors && right_neighbors                   %((c_temp(c_row-1,c_col)+c_temp(c_row-2,c_col)+c_temp(c_row-3,c_col)+c_temp(c_row-4,c_col))>pix_value)&&((c_temp(c_row,c_col+1)+c_temp(c_row,c_col+2)+c_temp(c_row,c_col+3)+c_temp(c_row,c_col+4))>pix_value)
                    c(c_row,c_col,:) = pix_value;

                elseif lower_neighbors && upper_neighbors                   %((c_temp(c_row+1,c_col)+c_temp(c_row+2,c_col)+c_temp(c_row+3,c_col)+c_temp(c_row+4,c_col))>pix_value)&&((c_temp(c_row-1,c_col)+c_temp(c_row-2,c_col)+c_temp(c_row-3,c_col)+c_temp(c_row-4,c_col))>pix_value)
                    c(c_row,c_col,:) = pix_value;

                elseif left_neighbors && right_neighbors                    %((c_temp(c_row,c_col-1)+c_temp(c_row,c_col-2)+c_temp(c_row,c_col-3)+c_temp(c_row,c_col-4))>pix_value)&&((c_temp(c_row,c_col+1)+c_temp(c_row,c_col+2)+c_temp(c_row,c_col+3)+c_temp(c_row,c_col+4))>pix_value)
                    c(c_row,c_col,:) = pix_value;

                elseif left_neighbors && lower_neighbors                    %((c_temp(c_row,c_col-1)+c_temp(c_row,c_col-2)+c_temp(c_row,c_col-3)+c_temp(c_row,c_col-4))>pix_value)&&((c_temp(c_row+1,c_col)+c_temp(c_row+2,c_col)+c_temp(c_row+3,c_col)+c_temp(c_row+4,c_col))>pix_value)
                    c(c_row,c_col,:) = pix_value;

                elseif right_neighbors && lower_neighbors                   %((c_temp(c_row,c_col+1)+c_temp(c_row,c_col+2)+c_temp(c_row,c_col+3)+c_temp(c_row,c_col+4))>pix_value)&&((c_temp(c_row+1,c_col)+c_temp(c_row+2,c_col)+c_temp(c_row+3,c_col)+c_temp(c_row+4,c_col))>pix_value)
                    c(c_row,c_col,:) = pix_value;

                %else
                %    c(c_row,c_col,:) = 0.0;                                 %Is this else unnecessary? If we get into this loop in the first place, since c was initialized as zeros, 
                end
            end
        end
    end
    c_temp=c(:,:,1);
end

%If using new geometry, generate two vector fields: pixel positions
%(relative to the center of Mars), and unit surface normals. There are
%different methods to do this and we are still trying to find the best one.
if geoflags.offnadir
   pix_positions = zeros(c_num_rows,c_num_cols,3);
   pix_normals = zeros(c_num_rows,c_num_cols,3);
   if geoflags.flat
       for c_row = 1:c_num_rows
           theta = pi/2-ddr_lat(c_row);
           for c_col = 1:c_num_cols
               phi = ddr_lon(c_col);
               pix_normals(c_row,c_col,:) = reshape([sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)], [1,1,3]);
               %For the flat approximation, we'll say all pixels have equal
               %elevations, equal to the mean of the DDR elevations
               pix_positions(c_row,c_col,:) = (scene_areoid+mean(mean(ddr_iof(:,:,10))))*reshape(pix_normals(c_row,c_col,:),[1,1,3]);
           end
       end 
   else
        %Below are various interpolation functions we have tried.
        %(1) Use the MATLAB function griddata() 
        %[pix_positions, pix_normals] = normal_spline(ddr_iof,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_spacing,scene_areoid,c(:,:,1));
        %(2) Interpolate with a homemade Gaussian filter. Really jagged so far.
        %[pix_positions, pix_normals] = gaussian_interp(ddr_iof,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_spacing,scene_areoid,c(:,:,1));
        %(3) Don't interpolate, and just assume it's the same as the target
        %pixel during kernel calculation.
        %[pix_positions, pix_normals] = nointerp(ddr_iof,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_spacing,scene_areoid,c(:,:,1));
        %(4) If we happen to have HiRISE coverage of the area, use that
        %high-resolution grid.
        %[pix_positions, pix_normals] = hirise_read(ddr_iof,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_spacing,scene_areoid,c(:,:,1),filenames.hirise_filename);
        %(5) Use the MATLAB function scatteredInterpolant
        [pix_positions, pix_normals] = triscat_interp(ddr_iof,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_spacing,scene_areoid,c(:,:,1));
   end
else
    %Still must return these variables, but we won't use them so we'll take
    %up as little memory as possible
    pix_normals = 0;
    pix_positions = 0;
end