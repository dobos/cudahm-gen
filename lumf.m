function phi = lumf(beta, l, u, x)
%--------------------------------------------------------------
% Calculates the value of the luminosity function
%--------------------------------------------------------------

% Low cut-off version
phi = (x ./ u) .^ beta .* (1 - exp(-x ./ l)) .* exp(-x ./ u);

% Log version
%phi = beta .* log(x ./ u) - x ./ u + log(1 - exp(-x ./ l));

% Original Schechter form
%phi = (x ./ u) .^ beta .* exp(-x ./ u);