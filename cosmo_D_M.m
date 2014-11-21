function D_M = cosmo_D_M(z, cosmo)
%--------------------------------------------------------------
% Calculates the comoving distance (transverse)
%--------------------------------------------------------------

D_C = cosmo_D_C(z, cosmo);

if cosmo.Omega_K > 0
    % hyperbolic
    D_M = cosmo.D_Hubble * 1 / sqrt(cosmo.Omega_K) * sinh(sqrt(cosmo.Omega_K) * D_C / cosmo.D_H);
elseif cosmo.Omega_K == 0
    % flat
    D_M = D_C;
elseif cosmo.Omega_K < 0
    % spherical
    D_M = cosmo.D_Hubble * 1 / sqrt(-cosmo.Omega_K) * sin(sqrt(-cosmo.Omega_K) * D_C / cosmo.D_H);
end