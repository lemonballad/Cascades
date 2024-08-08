function prefactor=prefactor_3_5_offres(l,C_Solv,alpha2,mu_eg,n_w_t_Solv,w_t)
% Computes the concentration dependent prefactor for off-resonant resonant
% mixed experiment
%
%   l                  : Path length in m
%   C_Solv        : Concentraion of solvent in mol/m^3
%   alpha2        : Polarizability squared of solvent in  A^4 s^8 kg^-2
%   mu_eg        : Magnitude transition dipole in D
%   n_w_t_Sol  : Refractive index at signal frequency w_t
%   w_t             : signal frequency in wavenumbers

%% Define constants
c=2.998E-7; % speed of light in m/fs
D_SI=3.3356E-30; % Conversion factor from Debye to SI units: D => s A m
e_0=8.8542E-12; % electric permitivity in a vaccum in s^4 A^2 m^-3 kg^-1
hbar=1.0546E-34; % reduced Planck constant in J s => kg m^2 s^-1
m_cm=1E2; % Conversion factor from meters to centimeters cm m^-1
N_A=6.022E23; % Avogadro's number in mol^-1

%% Do necessary conversions
% Convert mu_eg to SI units
mu_eg=mu_eg*D_SI; % (D) (s A m D^-1) => s A m

% Convert concentration to number density C => N
N_Solv=C_Solv*N_A; % (mol m^-3) (mol^-1) => m^-3

% Convert signal frequency to fHz^-1
w_t=w_t*c*2*pi*m_cm;

%% Compute the terms of the prefactor
% Orientational constant
prefactor=7/25;

% Second term
prefactor=prefactor*(l*N_Solv*w_t)/(n_w_t_Solv*c); % ((m) (m^-3) (fs^-1)) (m fs^-1)^-1 => m^-3

% Third term
prefactor=prefactor*hbar/(2*e_0)*alpha2/mu_eg^2; % (m^-3) (kg m^2 s^-1)  (A^4 s^8 kg^-2) (s A m)^-2 (s^4 A^2 m^-3 kg^-1)^-1 => s

end % End function prefactor_3_5