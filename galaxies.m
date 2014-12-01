function res = galaxies(n, maglim, rz, rl)
%--------------------------------------------------------------
% Generates a random sample of galaxies
%--------------------------------------------------------------

sz = size(rz.z);
nz = sz(1);

sl = size(rl);
nl = sl(1);

% Preallocate results
res = struct( ...
    'L', zeros([n, 1]), ...
    'z', zeros([n, 1]), ...
    'D_L', zeros([n, 1]), ...
    'F', zeros([n, 1]), ...
    'F_err', zeros([n, 1]), ...
    'F_real', zeros([n, 1]), ...
    'AbsM', zeros([n, 1]), ...
    'm', zeros([n, 1]), ...
    'm_err', zeros([n, 1]), ...
    'm_real', zeros([n, 1]));

for i = 1:n
    
    while true
        % Pick a redshift from the pre-generated sample
        zi = int32(floor(1 + nz * rand()));
        z = rz.z(zi);

        % Pick a luminosity from the pre-generated sample
        li = int32(floor(1 + nl * rand()));
        lm = rl(li);

        % Calculate the absolute and apparent magnitude
        absmag = -2.5 * log10(lm);
        D_L = rz.D_L(zi);
        DM = 5 * log10(D_L) + 25;
        mag = absmag + DM;
        
        % Use realistic error based on SDSS photometry
        pa = [0.00000003, 0, 0.000001, 0, 0.00005, 0, 0.0015];
        err = polyval(pa, mag - 12);
        pa = [0.000003, 0, 0, 0, 0];
        err = err + abs(normrnd(0, polyval(pa, mag - 12)));
        
        % Add error
        magerr = mag + err * sign(rand() - 0.5);

        % Only accept galaxy within flux limit, otherwise
        % sample from luminosity function
        if (magerr > maglim)
            continue
        end

        % Now we have a valid pair of z and luminosity

        % Convert to flux
        fl = 10 .^ (-0.4 * magerr);
        fl_real = 10 .^ (-0.4 * mag);
        fl_err = abs(fl - fl_real);
        
        % Copy into results
        res.L(i) = lm;
        res.z(i) = z;
        res.D_L(i) = D_L;
        res.F(i) = fl;
        res.F_err(i) = fl_err;
        res.F_real(i) = fl_real;
        res.AbsM(i) = absmag;
        res.m(i) = magerr;
        res.m_err(i) = err;
        res.m_real(i) = mag;        
        break;
    end
end
