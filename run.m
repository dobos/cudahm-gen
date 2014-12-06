% Luminosities are pre-generated for faster debugging

cosmo = cosmo_init();   % init standard cosmology

% Sampled volume parameters
z_min = 0.01;           % minimum redshift
z_max = 0.1;            % maximum redshift

m_limit = 18;           % Flux limit in apparent magnitude

%mmin = -15;
mmin = m_limit - cosmo_DM(z_min, cosmo);
mmax = -23;

% Generate luminosities
N = 10000;
beta = -2.2;
lmag = -16;
umag = -20;
rl = luminosities(N, beta, lmag, umag, mmin, mmax);

% Read random luminosities generated previously
% rl = dlmread('lm_betalu.dat');

% Intergate luminosity function
I = integral(@(x)lumf(beta, mag2lum(lmag), mag2lum(umag), x), mag2lum(mmin), mag2lum(mmax));

% Plot generated sample and luminosity function

    % lin-lin plot (not much to see)
    figure;
    [n xout] = hist(rl, 100);
    bar(xout, n, 'barwidth', 1, 'basevalue', 0); 
    hold;
    plot(xout, lumf(beta, mag2lum(lmag), mag2lum(umag), xout), 'r');
    xlabel('L'); ylabel('phi');

    % lin-log plot
    figure;
    [n xout] = hist(rl, 100);
    bar(xout, n, 'barwidth', 1, 'basevalue', 0); 
    hold;
    plot(xout, lumf(beta, mag2lum(lmag), mag2lum(umag), xout), 'r');
    set(gca,'YScale','log'); xlabel('L'); ylabel('phi');

    % log-log plot (mags on x axis)
    figure;
    [n xout] = hist(lum2mag(rl), 50);
    bar(xout, log10(n), 'barwidth', 1, 'basevalue', 0);
    hold;
    plot(xout, 1 + log10(lumf(beta + 1, mag2lum(lmag), mag2lum(umag), mag2lum(xout))), 'r');
    xlabel('M'); ylabel('log phi');

% Generate redshifts
rz = redshifts(50000, z_min, z_max);

gals = galaxies(10000, m_limit, rz, rl);

% Plot galaxy magnitude distribution

    figure;
    [n xout] = hist(gals.AbsM, 50);
    bar(xout, log10(n), 'barwidth', 1, 'basevalue', 0);
    hold;
    plot(xout, 2.5 + log10(lumf(beta + 1, mag2lum(lmag), mag2lum(umag), mag2lum(xout))), 'r');
    xlabel('M'); ylabel('log phi');
        
% Calculate absolute magnitude (luminosity) limit from apparent magnitude
% (flux) limit as a function of redshift

zz = z_min:0.0001:z_max;
AbsM_limit = m_limit - cosmo_DM(zz, cosmo);

% Plot absolute magnitudes and limiting function
    figure;
    scatter(gals.z, gals.AbsM, 1, '.');
    hold;
    plot(zz, AbsM_limit, 'r');
    
% Weighting function comes from the V_max method
% For every AbsMag, calculate the maximum redshift at which it falls
% out of the magnitude limit, this is the inverse of AbsM_limit(z) as
% computed above.

% Plot the maximum redshift visible as a function of absolute magnitude
    figure;
    plot(AbsM_limit, zz);

% Now convert z to volume

V_min = cosmo_V_C(z_min, cosmo);
V_max = cosmo_V_C(z_max, cosmo);

V_limit = cosmo_V_C(zz, cosmo) - V_min;

% Plot the volume visible as a function of absolute magnitude
    figure;
    plot(AbsM_limit, V_limit);
    
% Compute the weighting function, normalize for 1
% To fit W as a function of L, we only calculate it in the meaningful
% range where the flux limit applies. Outside this range the value is 1

mm = m_limit - cosmo_DM(z_max, cosmo):0.1:m_limit - cosmo_DM(z_min, cosmo);
W = interp1(AbsM_limit, V_limit ./ (V_max - V_min), mm);

% Zero out abs mags that are out of the limit even at z_min
% this part of the function won't be used
% W(mm > m_limit - cosmo_DM(z_min, cosmo)) = 0;

% Set weight to one if galaxy with abs mag is visible everywhere
% W(mm < m_limit - cosmo_DM(z_max, cosmo)) = 1;

% Plot W as a function of abs mag
    figure;
    plot(mm, W);
    
% Fit W with a polynom in the interesting luminosity range
p = polyfit(mag2lum(mm), W, 2));

    % Plot W as a function of L
    figure;
    plot(mag2lum(mm), W);
    hold;
    plot(mag2lum(mm), polyval(p, mag2lum(mm)), 'r');