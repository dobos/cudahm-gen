function D_C = cosmo_D_C(z, cosmo)
%--------------------------------------------------------------
% Calculates the comoving distance
%--------------------------------------------------------------

% Comoving distance
D_C = cosmo.D_Hubble * arrayfun(@(y)integral(@(x)cosmo_E(x,cosmo), 0, y), z);