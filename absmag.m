function M = absmag(mag, z, cosmo)
%--------------------------------------------------------------
% Converts apparent magnitude to absolute
%--------------------------------------------------------------

M = mag - DM(z, cosmo);