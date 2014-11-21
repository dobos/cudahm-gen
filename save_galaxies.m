function save_galaxies(galaxies, filename)
%--------------------------------------------------------------
% Writes a sample of galaxies into files
%--------------------------------------------------------------

sz = size(galaxies.m);
n = sz(1);

ff = strcat('filtered_fluxes_', filename, '.dat');
dlmwrite(ff, [galaxies.F zeros([n 1])], ' ');

ff = strcat('fluxes_', filename, '.dat');
dlmwrite(ff, galaxies.F_real, ' ');

ff = strcat('dist_', filename, '.dat');
dlmwrite(ff, galaxies.D_L, ' ');