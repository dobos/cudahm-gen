function r = randdist(func, a, b)
%--------------------------------------------------------------
% Generates a random number between a and b from the
% distribution @func. @func must be positive and monothonic.
%--------------------------------------------------------------

m = max([func(a) func(b)]);

% make sure first iteration executes
px = 0;
y = 1;

while y > px
    r = a + (b - a) * rand();
    y = m * rand();
    px = func(r);
end