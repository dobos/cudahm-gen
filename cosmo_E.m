function E = cosmo_E(z, cosmo)
%--------------------------------------------------------------
% Calculates the value of E(z) for cosmological distances
%--------------------------------------------------------------

E = cosmo.Omega_R * (1 + z) .^ 4 + ...
    cosmo.Omega_M * (1 + z) .^ 3 + ...
    cosmo.Omega_K * (1 + z) .^ 2 + cosmo.Omega_L;

E = 1 ./ sqrt(E);