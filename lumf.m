function phi = lumf(beta, l, u, x)
%--------------------------------------------------------------
% Calculates the value of the luminosity function
%--------------------------------------------------------------

phi = (x ./ u) .^ beta .* (1 - exp(-x ./ l)) .* exp(-x ./ u);
%phi = (x ./ u) .^ beta .* exp(-x ./ u);