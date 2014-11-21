function lum = mag2lum(mag)
%--------------------------------------------------------------
% Converts magnitude to luminosity
%--------------------------------------------------------------

lum = 10 .^ (-0.4 * mag);