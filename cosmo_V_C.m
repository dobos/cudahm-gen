function V_C = cosmo_V_C(z, cosmo)
%--------------------------------------------------------------
% Calculates the comoving volume of the universe
%--------------------------------------------------------------

D_M = cosmo_D_M(z, cosmo);

if cosmo.Omega_K > 0
    % hyperbolic
    V_C = -1;   % TBW
elseif cosmo.Omega_K == 0
    % flat
    V_C = 4 * pi / 3 * D_M .^ 3;
elseif cosmo.Omega_K < 0
    % spherical
    V_C = -1;   % TBW
end