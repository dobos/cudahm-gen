function lm = luminosities(n, beta, lmag, umag, m1, m2)
%--------------------------------------------------------------
% Generates a random sample based on the luminosity function
%--------------------------------------------------------------

l = mag2lum(lmag);
u = mag2lum(umag);
l1 = mag2lum(m1);
l2 = mag2lum(m2);

% Convert mags to luminosity
plmax = lumf(beta, l, u, mag2lum(m1));   % max at faint end

% Preallocate
lm = zeros([n, 1]);

datestr(now, 'HH:MM:SS')

for i = 1:n

    y = 1;
    p = 0;
    
    while y > p
        % Pick a luminosity between l1 and l2
        rl = l1 + (l2 - l1) * rand();
        
        % Pick a magnitude between m1 and m2
        %rm = m1 + (m2 - m1) * rand();
        %rl = mag2lum(rm);
        
        p = lumf(beta, l, u, rl);
        y = plmax * rand();
    end
    
    lm(i) = rl;
end

datestr(now, 'HH:MM:SS')