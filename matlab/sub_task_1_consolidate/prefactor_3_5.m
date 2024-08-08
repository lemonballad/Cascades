function prefactor=prefactor_3_5(l,C,mu_eg,n_w_t,w_t)
% Computes the concentration dependent prefactor for eq. 22
% B. P. Molesky, P. G. Giokas, Z. Guo, A. M. Moran, "Multidimensional
% Resonance Raman Spectroscopy by Six-Wave Mixing in the Deep UV" J. Chem.
% Phys., 141, 114202 (2014).
%
%   l               : Path length in m
%   C              : Concentraion in mol/L
%   mu_eg     : Magnitude transition dipole in D
%   n_w_t      : Refractive index at signal frequency w_t
%   w_t          : signal frequency in wavenumbers

%% Define constants
c=2.998E-7; % speed of light in m/fs
D_SI=3.3356E-30; % Conversion factor from Debye to SI units: D => s A m
e_0=8.8542E-12; % electric permitivity in a vaccum in s^4 A^2 m^-3 kg^-1
hbar=1.0546E-34; % reduced Planck constant in J s => kg m^2 s^-2
L_m3=1E3; % Conversion factor from liters to meters cubed
m_cm=1E2; % Conversion factor from meters to centimeters
N_A=6.022E23; % Avogadro's number in mol^-1

%% Do necessary conversions
% Convert mu_eg to SI units
mu_eg=mu_eg*D_SI; % (D) (s A m D^-1) => s A m

% Convert concentration to number density C => N
N=C*N_A*L_m3; % (mol L^-1) (mol^-1) (L m^-3) => m^-3

% Convert signal frequency to fHz^-1
w_t=w_t*c*2*pi*m_cm;

%% Compute the terms of the prefactor
% Orientational constant
prefactor=7/25;

% Second term
prefactor=prefactor*(l*N*w_t)/(n_w_t*c); % ((m) (m^-3) (fs^-1)) (m fs^-1)^-1 => m^-3

% Third term
prefactor=prefactor*(mu_eg^2)/(2*e_0*hbar); % (m^-3) (s A m) ((s^4 A^2 m^-3 kg^-1) (kg m^2 s^-2))^-1 => s^-1

end % End function prefactor_3_5