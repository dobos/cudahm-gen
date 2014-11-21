function z = redshifts(n, z1, z2)
%--------------------------------------------------------------
% Generates a random sample of redshifts, assuming uniform
% galaxy density in the universe
%--------------------------------------------------------------

cosmo = cosmo_init();
pzmax = cosmo_dV_C(z2, cosmo);  % max at high z

% Preallocate
z = struct( ...
    'z', zeros([n, 1]), ...
    'D_L', zeros([n, 1]));

for i = 1:n
    % Pick a redshift between z1 and z2
    y = 1;
    p = 0;
    
    while y > p
        rz = z1 + (z2 - z1) * rand();
        cc = cosmo_all(rz, cosmo);
        p = cc.dV_C;
        y = pzmax * rand();
    end
    
    z.z(i) = rz;
    z.D_L(i) = cc.D_L;
end
