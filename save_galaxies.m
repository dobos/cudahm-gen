function save_galaxies(galaxies, path)
%--------------------------------------------------------------
% Writes a sample of galaxies into files
%--------------------------------------------------------------

sz = size(galaxies.m);
n = sz(1);

ff = strcat(path, '/filtered_fluxes.dat');
dlmwrite(ff, [galaxies.F galaxies.F_err], ' ');

ff = strcat(path, '/fluxes.dat');
dlmwrite(ff, galaxies.F_real, ' ');

ff = strcat(path, '/dist.dat');
dlmwrite(ff, galaxies.D_L, ' ');