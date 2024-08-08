function [inhomog,gamma_0,kappa,lambda,weg,translen,n_w_t_Sol,...
    wvib,disp,K,rho_Solv,mm_Solv,wsolv,n_w_t_Solv,weg_Solv]...
    =DEFENSE_PNA_parameters(solvent)
% Returns spectral simulation parameters of p-nitroaniline
%   solvent             : the solvent parameters requested for
%
% Solvent options include: cyclohexane, 1,4-dioxane, dichloromethane,
% acetonitrile, and methanol.
%
% If one of the above options is not used, methanol is chosen
%
% A. M. Moran, A. M. Kelley. "Solvent effects on ground and excited
% electronic state structures of p-nitroaniline."

% Determine solvent
switch solvent
    case 'cyclohexane'
        % Electronic inhomog std dev. (cm^-1)
        inhomog=0;
        % Electronic homog (FWHM) (cm^-1)
        gamma_0=2045;
        % LAMBDA/DELTA
        kappa=0.1;
        % Classic reorganization energy: DEL^2/(2kT) (cm^-1)
        lambda=1890;
        % Electronic origin (cm^-1
        weg=27542;
        % transition length (Angs)
        translen=1.1;
        % Refractive index
        n_w_t_Sol=1.427;
        % Vibrational frequencies (cm^-1)
        wvib=[859 1112 1338 1498 1599];
        % Dimensionless displacements
        disp=[0.88 0.288 1.433 0.382 0.496];
        % Solvent Raman differential cross section prefactor in m^2
        K=9.99e-31; 
        % Solvent density in g/m^3
        rho_Solv=778000; 
        % Solvent molar mass in g/mol
        mm_Solv=84.16; 
        % Solvent vibrational mode in cm^-1
        wsolv=801; 
        % refractive index Solvent
        n_w_t_Solv=1.43; 
        % Electronic energy gap origin Solvent in cm^-1
        weg_Solv=115000; 
    case '1,4-dioxane'
        inhomog=110;
        gamma_0=2473;
        kappa=0.05;
        lambda=2700;
        weg=24100;
        translen=1.189;
        n_w_t_Sol=1.422;
        wvib=[861 1113 1328 1508 1602];
        disp=[1.23 0.53 1.17 0.223 0.341];
        K=4.41e-31; 
        rho_Solv=1033000; 
        mm_Solv=88.11; 
        wsolv=835; 
        n_w_t_Solv=1.42; 
        weg_Solv=99600; 
    case 'dichloromethane'
        inhomog=60;
        gamma_0=2915;
        kappa=0.07;
        lambda=3780;
        weg=23490;
        translen=1.172;
        n_w_t_Sol=1.424;
        wvib=[862 1112 1326 1504 1602 1625];
        disp=[1.049 0.231 1.184 0.1727 0.3195 0.27];
        K=1.92e-30; 
        rho_Solv=1327000; 
        mm_Solv=84.93; 
        wsolv=702; 
        n_w_t_Solv=1.42; 
        weg_Solv=126900; 
    case 'acetonitrile'
        inhomog=0;
        gamma_0=3000;
        kappa=0.2;
        lambda=4310;
        weg=22720;
        translen=1.199;
        n_w_t_Sol=1.344;
        wvib=[861 1112 1326 1507 1599];
        disp=[0.84 0.242 0.976 0.152 0.292];
        K=8.99e-32; 
        rho_Solv=786000; 
        mm_Solv=41.05; 
        wsolv=919; 
        n_w_t_Solv=1.34; 
        weg_Solv=105700; 
    case 'methanol'
        inhomog=750;
        gamma_0=1410;
        kappa=0.1;
        lambda=900;
        weg=23800;
        translen=1.243;
        n_w_t_Sol=1.329;
        wvib=[861 1113 1328 1508 1600];
        disp=[1.46 0.65 1.451 0.266 0.287];
        K=7.62e-31; 
        rho_Solv=792000; 
        mm_Solv=32.04; 
        wsolv=1035; 
        n_w_t_Solv=1.34; 
        weg_Solv=153100; 
    otherwise
        error('Error: This is not a valid solvent');
end
end % End function PNA_parameters