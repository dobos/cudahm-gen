function distmod = cosmo_DM(z, cosmo)
%--------------------------------------------------------------
% Calculates the distance modulus
%--------------------------------------------------------------

distmod = 5 * log10(cosmo_D_L(z, cosmo)) + 25;