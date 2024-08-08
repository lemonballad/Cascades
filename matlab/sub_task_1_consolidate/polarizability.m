function alpha2=polarizability(K,weg_Solv,w_L)
% Compute polarizability squared from Raman differential cross section
%   K                       : Raman differential cross section prefactor in
%                              m^2
%   weg_Solv         : First solvent electronic energy gap in cm^-1
%   w_L                  : Laser excitation frequency in cm^-1

%% Define constants
% Permitivity of free space in A^2 s^4 kg^-1 m^-3
e0=8.85418782e-12;
% Speed of light in m s^-1
cms=2.99792458e8;
% Speed of light in cm s^-1
ccms=cms*100;

%% Change units
% Convert units from cm^-1 => s^-1
weg_Solv=weg_Solv*ccms;
% Convert units from cm^-1 => s^-1
w_L=w_L*ccms;

%%
% Compute the frequency dependence of the polarizability in s^4
wfactor=((weg_Solv^2+w_L^2)/(weg_Solv^2-w_L^2)^2)^2;
% Multiply constants in A^4 s^4 kg^-2 m^-2
pfactor=16*pi^2*e0^2*cms^4;
% Polarizability squared in A^4 s^8 kg^-2
alpha2=pfactor*K*wfactor;

return % End function polarizability
