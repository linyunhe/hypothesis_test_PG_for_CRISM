function [ projected_scene, mapProjectionInfo ] = crism_map_project( c, projection, ddr_iof, ddr_lat, ddr_lon, pix_size )
%CRISM_MAP_PROJECT Produce Mars Equirectangular map from C.
%   Converts the local-Cartesian image C to an official map projeciton.
%   The output has North upward (low-numbered rows) and East to the right
%   (high-numbered cols), unlike c.
%
%   ddr_lat and ddr_lon contain lat & lon info for c.
%
%   Supported projections are currently:
%   - 'Mars_Equirect': Equirectangular.
%   - 'None': None. Just returns Local-Cartesian scene flipped around so
%       that it's suitable for saving.
%

r = 3396190.0; %radius of Mars in meters (exact value used in SuperGLT reconstruction)
[c_num_rows, c_num_cols, num_bands] = size(c);

% Project the Scene (different ways depending on requested projection)
if strcmp(projection,'Mars_Equirect')
    % Mars Equirectangular: pixels with uniform delta-latitudes &
    % delta-longitudes. Definition:
    % X = R * longitude
    % Y = R * latitude
    %
    % We set scene up so our reference pixel (highest lat, lowest lon) has
    % size pix_size * pix_size m^2.
    % For the reference pixel, delta-lon = delta-X/R = pix_size/R
    % This becomes our constant scaling of latitudes and longitudes
    % throughout the projected c.
    
    min_lat = min(min(ddr_iof(:,:,4)*pi/180.0));
    max_lat = max(max(ddr_iof(:,:,4)*pi/180.0));
    min_lon = min(min(ddr_iof(:,:,5)*pi/180.0));
    max_lon = max(max(ddr_iof(:,:,5)*pi/180.0));
    
    delta_rad = pix_size / r; % The size in both directions of all pixels in projected scene
    
    proj_num_rows = 1+round((max_lat-min_lat) / delta_rad);
    proj_num_cols = 1+round((max_lon-min_lon) / delta_rad);
    
    projected_scene = zeros(proj_num_rows,proj_num_cols,num_bands);
    
    proj_lat_spacing = abs(max_lat-min_lat)/(proj_num_rows-1);
    proj_lon_spacing = abs(max_lon-min_lon)/(proj_num_cols-1);
    
    %determine latitude for every projected row given that projected rows are evenly spaced
    proj_lats_small = max_lat-proj_lat_spacing*((1:proj_num_rows)-1);

    %determine longitude for every projected col given that projected columns are evenely spaced
    proj_lons_small = min_lon+proj_lon_spacing*((1:proj_num_cols)-1);
    
    % preparing C to be interpolated
    c_lats = repmat(ddr_lat',[1 c_num_cols]);
    c_lons = repmat(ddr_lon,[c_num_rows 1]);
    
    proj_lats = repmat(proj_lats_small',[1 proj_num_cols]);
    proj_lons = repmat(proj_lons_small,[proj_num_rows 1]);
    
    % Interpolation. Form each pixel in the projected scene by
    % interpolating from C.
    for band = 1:num_bands
        projected_scene(:,:,band) = interp2(c_lons,c_lats,c(:,:,band),proj_lons,proj_lats);
    end

elseif strcmp(projection,'None')
    projected_scene = transformScene( c ); % double-flip for correct directionality
    
else
    disp('Unknown map projection requested!');
    projected_scene = [];
    
end

% Gather map projection info for header file
if strcmp(projection,'Mars_Equirect')
    equirectMPerDeg = r * pi/180; % using r assumed by projection

    % tie point
    % uses upper-left corner with c's lower-right corner's geocoordinates to
    % compensate for double flip applied to make transformed scene
    mapProjectionInfo.tieRow = 1+.5;
    mapProjectionInfo.tieLatM = ddr_lat(c_num_rows) * 180/pi * equirectMPerDeg;
    mapProjectionInfo.tieCol = 1+.5;
    mapProjectionInfo.tieLonM = ddr_lon(c_num_cols) * 180/pi * equirectMPerDeg;

    % spacing & rotation
    mapProjectionInfo.pixelSizeX = pix_size;
    mapProjectionInfo.pixelSizeY = pix_size;
    mapProjectionInfo.rotation = 0; % north is up
    
elseif strcmp(projection,'None')
    equirectMPerDeg = r * pi/180; % using r assumed by projection

    % tie point
    % uses upper-left corner with c's lower-right corner's geocoordinates to
    % compensate for double flip applied to make transformed scene
    mapProjectionInfo.tieRow = 1+.5;
    mapProjectionInfo.tieLatM = ddr_lat(c_num_rows) * 180/pi * equirectMPerDeg;
    mapProjectionInfo.tieCol = 1+.5;
    mapProjectionInfo.tieLonM = ddr_lon(c_num_cols) * 180/pi * equirectMPerDeg;

    % spacing & rotation
    mapProjectionInfo.pixelSizeX = pix_size;
    mapProjectionInfo.pixelSizeY = pix_size;
    mapProjectionInfo.rotation = 0; % north is up
    
else
    disp('Unknown map projection requested!');
    projected_scene = [];
end

end

