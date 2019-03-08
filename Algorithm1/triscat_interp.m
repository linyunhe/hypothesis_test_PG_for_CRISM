function [pix_positions,pix_normals] = triscat_interp(ddr_iof,c_num_rows,c_num_cols,ddr_lat,ddr_lon,lat_number,lon_number,pix_spacing,scene_areoid,interp_scene);

body499_radii = [3396.19,3396.19,3376.20]; %from the planetary constants kernel in SPICE 
size_kernel = 21;
padding = (size_kernel-1)/2;

[mu_num_rows, mu_num_cols, ~] = size(ddr_iof);
num_samples = mu_num_rows*mu_num_cols;

orig_radii = NaN*ones(c_num_rows,c_num_cols);
orig_normals = NaN*ones(c_num_rows,c_num_cols,3);

%The following loop is for comparison purposes only
for mu_row = 1:mu_num_rows
    for mu_col = 1:mu_num_cols
        
        c_row = lat_number(mu_row, mu_col);
        c_col = lon_number(mu_row, mu_col);
        
        orig_radii(c_row,c_col) = scene_areoid+ddr_iof(mu_row,mu_col,10);
        
        lat = ddr_iof(mu_row,mu_col,4)*pi/180;
        theta = pi/2-lat;
        lon = ddr_iof(mu_row,mu_col,5)*pi/180;
        phi = lon;
        r = scene_areoid + ddr_iof(mu_row, mu_col, 10);
        rhat = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
        gradient_factor = 2./(body499_radii.^2); 
        ellipsoid_normal = gradient_factor.*rhat;
        ellipsoid_normal = ellipsoid_normal/sqrt(ellipsoid_normal*ellipsoid_normal');
        %find surface normal vector
        %first, get the northward direction (compass pointing direction)
        zhat = [0,0,1];
        tangent_north = zhat-(zhat*ellipsoid_normal')*ellipsoid_normal;
        %find dip (azimuth) direction
        azimuth = ddr_iof(mu_row,mu_col,9)*pi/180;
        %dip direction is in the plane of tangent_north and
        %cross(tangent_north, ellipsoid_normal), so it is just a linear combination:
        dip_direction = tangent_north*cos(azimuth) + cross(tangent_north, ellipsoid_normal)*sin(azimuth);
        dip_direction = dip_direction/sqrt(dip_direction*dip_direction');%normalization

        
        %find surface normal
        slope_angle = ddr_iof(mu_row,mu_col,8)*pi/180;
        %surf_normal is in the plane of ellipsoid_normal and dip_direction, so it is
        %just a linear combination
        surf_normal = ellipsoid_normal*cos(slope_angle) + dip_direction*sin(slope_angle); 
        surf_normal = surf_normal/sqrt(surf_normal*surf_normal'); %normalization
        
        orig_normals(c_row,c_col,:) = reshape(surf_normal,[1,1,3]);
    end
end

x = zeros(num_samples,1); %x will be rows
y = zeros(num_samples,1); %y will be columns
z = zeros(num_samples,1); %z will be elevations

sample = 0;
for mu_row = 1:mu_num_rows
    for mu_col = 1:mu_num_cols
        sample = sample+1;
        x(sample) = lon_number(mu_row,mu_col);
        y(sample) = lat_number(mu_row,mu_col);
        z(sample) = ddr_iof(mu_row,mu_col,10); %elevation
    end
end

%invoke the scatteredinterpolant class
F = scatteredInterpolant(x,y,z);
[xq,yq] = meshgrid(1:c_num_cols,1:c_num_rows);
zq = F(xq,yq);

%Now project this data as a map
pix_radii = scene_areoid + zq;
pix_positions = zeros(c_num_rows,c_num_cols,3);
for c_row = 1:c_num_rows
    lat = ddr_lat(c_row);
    theta = pi/2-lat;
    for c_col = 1:c_num_cols
        lon = ddr_lon(c_col);
        phi = lon;
        rhat = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
        pix_positions(c_row,c_col,:) = pix_radii(c_row,c_col)*reshape(rhat,[1,1,3]);
    end
end
%Finally, derive the surface normals
pix_normals = zeros(c_num_rows,c_num_cols,3);
for c_row = 1:c_num_rows
    for c_col = 1:c_num_cols
       minrow = max(c_row-1,1);
       mincol = max(c_col-1,1);
       maxrow = min(c_row+1,c_num_rows);
       maxcol = min(c_col+1,c_num_cols);
       lon_deriv = reshape(pix_positions(c_row,maxcol,:)-pix_positions(c_row,mincol,:),[1,3]);
       lat_deriv = reshape(pix_positions(maxrow,c_col,:)-pix_positions(minrow,c_col,:),[1,3]);
       %normalize
       lon_deriv = lon_deriv/sqrt(lon_deriv*lon_deriv');
       lat_deriv = lat_deriv/sqrt(lat_deriv*lat_deriv');
       normal = cross(lat_deriv,lon_deriv);
       normal = normal/sqrt(normal*normal');
       pix_normals(c_row,c_col,:) = reshape(normal,[1,1,3]);
    end
end

%Try an averaging filter
kernel = ones(size_kernel,size_kernel);

smoothed_normals = zeros(c_num_rows,c_num_cols,3);
for c_row = 1:c_num_rows
    for c_col = 1:c_num_cols
        kernel = zeros(size_kernel,size_kernel);
        minrow = max(c_row-padding,1);
        mincol = max(c_col-padding,1);
        maxrow = min(c_row+padding,c_num_rows);
        maxcol = min(c_col+padding,c_num_cols);
        
        k_minrow = padding+1-(c_row-minrow);
        k_mincol = padding+1-(c_col-mincol);
        k_maxrow = padding+1+(maxrow-c_row);
        k_maxcol = padding+1+(maxcol-c_col);
        
        for row_iter = minrow:maxrow
            k_row = k_minrow+(row_iter-minrow);
            for col_iter = mincol:maxcol
                k_col = k_mincol+(col_iter-mincol);
                dist_sq = (c_row-row_iter)^2+(c_col-col_iter)^2;
                kernel(k_row,k_col) = exp(-dist_sq*log(16)/(padding/2)^2);
            end
        end
        
        kernel_sum = sum(sum(kernel(k_minrow:k_maxrow,k_mincol:k_maxcol)));
        local_kernel = repmat(kernel/kernel_sum,[1,1,3]);
        smoothed_normals(c_row,c_col,:) = sum(sum(local_kernel(k_minrow:k_maxrow,k_mincol:k_maxcol,:).*pix_normals(minrow:maxrow,mincol:maxcol,:)));
    end
end

pix_normals = smoothed_normals;


end