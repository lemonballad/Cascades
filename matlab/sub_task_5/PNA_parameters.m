function [inhomog,gamma_0,kappa,lambda,weg,translen,refrac_index,...
    wvib,disp]=PNA_parameters(solvent)
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
        refrac_index=1.427;
        % Vibrational frequencies (cm^-1)
        wvib=[859 1112 1338 1498 1599];
        % Dimensionless displacements
        disp=[0.88 0.288 1.433 0.382 0.496];
    case '1,4-dioxane'
        inhomog=110;
        gamma_0=2473;
        kappa=0.05;
        lambda=2700;
        weg=24100;
        translen=1.189;
        refrac_index=1.422;
        wvib=[861 1113 1328 1508 1602];
        disp=[1.23 0.53 1.17 0.223 0.341];
    case 'dichloromethane'
        inhomog=60;
        gamma_0=2915;
        kappa=0.07;
        lambda=3780;
        weg=23490;
        translen=1.172;
        refrac_index=1.424;
        wvib=[862 1112 1326 1504 1602 1625];
        disp=[1.049 0.231 1.184 0.1727 0.3195 0.27];
    case 'acetonitrile'
        inhomog=0;
        gamma_0=3000;
        kappa=0.2;
        lambda=4310;
        weg=22720;
        translen=1.199;
        refrac_index=1.344;
        wvib=[861 1112 1326 1507 1599];
        disp=[0.84 0.242 0.976 0.152 0.292];
    case 'methanol'
        inhomog=750;
        gamma_0=1410;
        kappa=0.1;
        lambda=900;
        weg=23800;
        translen=1.243;
        refrac_index=1.329;
        wvib=[861 1113 1328 1508 1600];
        disp=[1.46 0.65 1.451 0.266 0.287];
    otherwise
        error('Error: This is not a valid solvent');
end
end % End function PNA_parameters