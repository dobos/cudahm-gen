function lm = luminosities(beta, lmag, umag, mag1, mag2, n)
%--------------------------------------------------------------
% Generates a random sample based on the luminosity function
%--------------------------------------------------------------

% Convert mags to luminosity
l = mag2lum(lmag);
u = mag2lum(umag);
l1 = mag2lum(mag1);
l2 = mag2lum(mag2);
plmax = lumf(beta, l, u, l1);   % max at faint end

% Preallocate
lm = zeros([n, 1]);

datestr(now, 'HH:MM:SS')

for i = 1:n
    % Pick a luminosity between l1 and l2
    y = 1;
    p = 0;
    
    while y > p
        rl = l1 + (l2 - l1) * rand();
        p = lumf(beta, l, u, rl);
        y = plmax * rand();
    end
    
    lm(i) = rl;
end

datestr(now, 'HH:MM:SS')