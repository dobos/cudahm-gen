function all = cosmo_all(z, cosmo)
%--------------------------------------------------------------
% Calculates the important cosmological quantities
%--------------------------------------------------------------

% Comoving distance
E = cosmo_E(z, cosmo);
D_C = cosmo.D_Hubble * arrayfun(@(y)integral(@(x)cosmo_E(x,cosmo), 0, y), z);

% Comoving distance (transverse), volume and derivative of volume
if cosmo.Omega_K > 0
    % hyperbolic
    D_M = cosmo.D_Hubble * 1 / sqrt(cosmo.Omega_K) * sinh(sqrt(cosmo.Omega_K) * D_C / cosmo.D_H);
    V_C = -1;   % TBW
    dV_C = -1;  % TBW
elseif cosmo.Omega_K == 0
    % flat
    D_M = D_C;
    V_C = 4 * pi / 3 * D_M .^ 3;
    dV_C = cosmo.D_Hubble .* 4 * pi * D_C .^ 2 ./ E;
elseif cosmo.Omega_K < 0
    % spherical
    D_M = cosmo.D_Hubble * 1 / sqrt(-cosmo.Omega_K) * sin(sqrt(-cosmo.Omega_K) * D_C / cosmo.D_H);
    V_C = -1;   % TBW
    dV_C = -1;  % TBW
end

% Luminosity disntance
D_L = (1 + z) .* D_M;

% Angular diameter distance
D_A = D_M / (1 + z);

all = struct( ...
    'z', z, ...
    'E', E, ...
    'D_C', D_C, ...
    'D_M', D_M, ...
    'D_A', D_A, ...
    'D_L', D_L, ...
    'V_C', V_C, ...
    'dV_C', dV_C);
