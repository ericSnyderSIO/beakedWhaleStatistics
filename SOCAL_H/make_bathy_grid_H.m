%% make_bathy_grid_H.m
% this file will take the GMRT bathymetry data and format it to a .mat file

% load the data
[X, Y, Z] = grdread2('F:\bathymetry\siteHmed2.grd');
surf(X,Y,Z)
% s.EdgeColor = 'none'

contour(X,Y,Z,'black','showtext','on');

save('F:\bathymetry\siteHgridMedium.mat','X','Y','Z')