function D_L = cosmo_D_L(z, cosmo)
%--------------------------------------------------------------
% Calculates the luminosity distance
%--------------------------------------------------------------

D_L = (1 + z) .* cosmo_D_M(z, cosmo);