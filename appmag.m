function m = appmag(Mag, z, cosmo)
%--------------------------------------------------------------
% Converts absolute magnitude to apparent
%--------------------------------------------------------------

m = Mag + DM(z, cosmo);