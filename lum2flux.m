function flux = lum2flux(lum, dist)
%--------------------------------------------------------------
% Converts luminosity to flux
%--------------------------------------------------------------

flux = lum / (4 * pi * dist .^ 2);