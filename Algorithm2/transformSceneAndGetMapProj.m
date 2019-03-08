function [ transformedScene, mapProjectionInfo ] = transformSceneAndGetMapProj( c, pix_spacing, ddr_lat, ddr_lon )
%TRANSFORMSCENEANDGETMAPPROJ Rotate Scene & build map proj info for it
%   Detailed explanation goes here

[c_num_rows c_num_cols ~] = size(c);

transformedScene = flipdim(flipdim(c,1),2); % flip over both spatial axes

% Build structure of information to use to put map projection info in
% output header files

equirectMPerDeg = (3396190.0 * 2 * pi) / 360; % using r assumed by projection

% tie point
% uses upper-left corner with c's lower-right corner's geocoordinates to
% compensate for double flip applied to make transformed scene
mapProjectionInfo.tieRow = 1+.5;
mapProjectionInfo.tieLatM = ddr_lat(c_num_rows) * 180/pi * equirectMPerDeg;
mapProjectionInfo.tieCol = 1+.5;
mapProjectionInfo.tieLonM = ddr_lon(c_num_cols) * 180/pi * equirectMPerDeg;

% spacing & rotation
mapProjectionInfo.pixelSizeX = pix_spacing;
mapProjectionInfo.pixelSizeY = pix_spacing;
mapProjectionInfo.rotation = 0; % north is up

end

