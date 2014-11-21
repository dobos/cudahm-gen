function D_A = cosmo_D_A(z, cosmo)
%--------------------------------------------------------------
% Calculates the angular diameter distance
%--------------------------------------------------------------

D_A = cosmo_D_M(z, cosmo) / (1 + z);