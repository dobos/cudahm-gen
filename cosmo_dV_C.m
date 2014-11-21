function dV_C = cosmo_dV_C(z, cosmo)
%--------------------------------------------------------------
% Calculates the differential of the comoving volume wrt z
%--------------------------------------------------------------

% Only for flat universe!

E = cosmo_E(z, cosmo);
D_C = cosmo_D_C(z, cosmo);
dV_C = cosmo.D_Hubble .* 4 * pi * D_C .^ 2 ./ E;