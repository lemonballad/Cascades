function [gamma_vib,gamma_eg,weg,wvib,disp,mu_eg,n_w_t,l,w_t]...
    =myoglobin_parameters
% Returns spectral simulation parameters of myoglobin
%   solvent             : the solvent parameters requested for
%
% B. P. Molesky, Z. Guo, T. P. Cheshire, A. M. Moran, "Two-Dimensional
% Resonance Raman Spectroscopy of Oxygen- and Water- Ligated Myoglobin" J.
% Chem. Phys., 145, .034203 (2016)

% Determine solvent
% Electronic inhomog std dev. (cm^-1)
gamma_vib=10;
% Electronic homog (FWHM) (cm^-1)
gamma_eg=750;
% Electronic origin (cm^-1
weg=23250;
% Vibrational frequencies (cm^-1)
wvib=[220 370 674 1356];
% Dimensionless displacements
disp=[0.47 0.20 0.26 0.34];
% Transition dipole (D)
mu_eg=8.8;
% Refractive index
n_w_t=1.39;
% Path length l (m)
l=0.22E-3;
% signal frequency (cm^-1)
w_t=weg;

end % End function myoglobin_parameters