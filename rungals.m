function rungals(path, N, z_min, z_max, m_limit, M_max, beta, lmag, umag)

    % Typical input parameters
    % z_min = 0.01;           % minimum redshift
    % z_max = 0.1;            % maximum redshift
    % m_limit = 18;           % Flux limit in apparent magnitude
    % M_max = -23;            % Absolute magnitude limit
    % N = 100000;             % Number of galaxies
    % beta = -2.2;
    % lmag = -16;
    % umag = -20;

    % Init standard cosmology
    cosmo = cosmo_init();   
    
    % Create folder for results
    mkdir(path);
    
    % Evaluate absolute magnitude limit at z_min, this will set the
    % lower limit of luminosity
    M_min = m_limit - cosmo_DM(z_min, cosmo);
    
    L_min = mag2lum(M_min);
    L_max = mag2lum(M_max);
    
    % Find magnitude and luminosity above which all galaxies are visible
    % in the entire volume
    M_cut = m_limit - cosmo_DM(z_max, cosmo);
    L_cut = mag2lum(M_cut);
    
    % Luminosities are pre-generated for faster debugging
    
        % Generate luminosities
        % rl = luminosities(100000, beta, lmag, umag, M_min, M_max);
        % dlmwrite('lmcache.dat', rl, ' ');

        % Read random luminosities generated previously
        rl = dlmread('lmcache.dat');
        [N_rl temp] = size(rl);

    % Intergate luminosity function
    I = integral(@(x)lumf(beta, mag2lum(lmag), mag2lum(umag), x), L_min, L_max);

    % Plot generated sample and luminosity function

        bins = 100;
        C = N_rl / bins / I * (L_max-L_min);
    
        % lin-lin plot (not much to see)
        img = figure;
        [n xout] = hist(rl, bins);
        bar(xout, n, 'barwidth', 1, 'basevalue', 0); 
        hold;
        plot(xout, C * lumf(beta, mag2lum(lmag), mag2lum(umag), xout), 'r');
        xlabel('L'); ylabel('phi');
        print(img, '-dpng', strcat(path, '/phi(L)_lin.png'));
        close(img);
        
        % lin-log plot
        img = figure;
        [n xout] = hist(rl, bins);
        bar(xout, n, 'barwidth', 1, 'basevalue', 0); 
        hold;
        plot(xout, C * lumf(beta, mag2lum(lmag), mag2lum(umag), xout), 'r');
        set(gca,'YScale','log'); xlabel('L'); ylabel('phi');
        print(img, '-dpng', strcat(path, '/log_phi(L).png'));
        close(img);

        % log-log plot (mags on x axis)
        img = figure;
        [n xout] = hist(lum2mag(rl), bins);
        bar(xout, log10(n), 'barwidth', 1, 'basevalue', 0);
        hold;
        plot(xout, log10(C) + log10(lumf(beta + 1, mag2lum(lmag), mag2lum(umag), mag2lum(xout))), 'r');
        xlabel('M'); ylabel('log phi');
        print(img, '-dpng', strcat(path, '/log_phi(M).png'));
        close(img);

    % Luminosities are pre-generated for faster debugging
        
        % Generate redshifts
        rz = redshifts(50000, z_min, z_max);
        dlmwrite('zcache.dat', rl, ' ');
        
        % Load redshifts from cache
        % rz = dlmread('zcache.dat');

    % Generate galaxies
        
    gals = galaxies(N, m_limit, rz, rl);

    % Plot galaxy magnitude distribution
        img = figure;
        [n xout] = hist(gals.AbsM, 50);
        bar(xout, log10(n), 'barwidth', 1, 'basevalue', 0);
        hold;
        plot(xout, 2.5 + log10(lumf(beta + 1, mag2lum(lmag), mag2lum(umag), mag2lum(xout))), 'r');
        xlabel('M'); ylabel('log phi');
        print(img, '-dpng', strcat(path, '/log_phi(M)_fluxlim.png'));
        close(img);

    % Calculate absolute magnitude (luminosity) limit from apparent magnitude
    % (flux) limit as a function of redshift

    zz = z_min:0.0001:z_max;
    AbsM_limit = m_limit - cosmo_DM(zz, cosmo);

    % Plot absolute magnitudes and limiting function
        img = figure;
        scatter(gals.z, gals.AbsM, 1, '.');
        hold;
        plot(zz, AbsM_limit, 'r');
        xlabel('z'); ylabel('M');
        print(img, '-dpng', strcat(path, '/M(z).png'));
        close(img);

    % Weighting function comes from the V_max method
    % For every AbsMag, calculate the maximum redshift at which it falls
    % out of the magnitude limit, this is the inverse of AbsM_limit(z) as
    % computed above.

    % Plot the maximum redshift visible as a function of absolute magnitude
        img = figure;
        plot(AbsM_limit, zz);
        xlabel('z'); ylabel('M_{limit}');
        print(img, '-dpng', strcat(path, '/M_limit(z).png'));
        close(img);

    % Now convert z to volume

    V_min = cosmo_V_C(z_min, cosmo);
    V_max = cosmo_V_C(z_max, cosmo);

    V_limit = cosmo_V_C(zz, cosmo) - V_min;

    % Plot the volume visible as a function of absolute magnitude
        img = figure;
        plot(AbsM_limit, V_limit);
        xlabel('M_{limit}'); ylabel('V_{limit}');
        print(img, '-dpng', strcat(path, '/V_limit(M_limit).png'));
        close(img);

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
        img = figure;
        plot(mm, W);
        xlabel('W'); ylabel('M');
        print(img, '-dpng', strcat(path, '/W(M).png'));
        close(img);

    % Fit W with a polynom in the interesting luminosity range
    p = polyfit(mag2lum(mm), W, 2);

        % Plot W as a function of L
        img = figure;
        plot(mag2lum(mm), W);
        hold;
        plot(mag2lum(mm), polyval(p, mag2lum(mm)), 'r');
        xlabel('W'); ylabel('L');
        print(img, '-dpng', strcat(path, '/W(L).png'));
        close(img);
        
    % Save data
    save_galaxies(gals, path);
    
    % Save parameters    
    ff = strcat(path, '/params.txt');
    fileID = fopen(ff,'w');

    fprintf(fileID,'z_min = %f\n', z_min);
    fprintf(fileID,'z_max = %f\n', z_max);
    fprintf(fileID,'m_limit = %f\n', m_limit);
    fprintf(fileID,'M_max = %f\n', M_max);
    fprintf(fileID,'L_min = %f\n', L_min);
    fprintf(fileID,'L_max = %f\n', L_max);
    fprintf(fileID,'L_cut = %f\n', L_cut);
    fprintf(fileID,'W = [ %e %e %e ]\n', p(1), p(2), p(3));

    fclose(fileID);