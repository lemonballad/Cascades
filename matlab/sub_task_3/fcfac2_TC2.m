function [ovlp] = fcfac2_TC(d,nquanta,wn,wm)
% Compute overlap integrals by Heller's Method
% Equation 7 Myers et al. JCP 77, 3857 (1982)
%
%   d                : unitless mode displacement
%   nquanta    : number of vibrational quanta
%   wn             : ground state vibrational frequency
%   wm            : excited state vibrational frequency

%% Calculate basic parameters
% Huang Rys factor
s=d^2/2;
% Ratio of ground and excited state frequencies
R=wn/wm;
% <0|0>
O_0_0=sqrt(2/(1+R))*R^0.25*exp(-s/(1+R));

%% 
for mm=1:nquanta
for nn=1:nquanta
m=mm-1; %actual quantum number for excited state
n=nn-1; %actual quantum number for ground state
%
num=0;
for k=0:min(mm,nn)-1;
O_mk_g=O_0_0 * ((R-1)/(R+1)).^(m/2-k/2) * hermiteH(m-k,-sqrt(2*s*R/(R^2-1)))...
    ./sqrt(2^(m-k)*factorial(m-k));

O_g_nk=O_0_0 * ((1-R)/(1+R)).^(n/2-k/2) * hermiteH(n-k,sqrt(2*s/(1-R^2)))...
    ./sqrt(2^(n-k)*factorial(n-k));

    num=num+(2*sqrt(R)/(1+R)).^k...
        .*1./factorial(k)...
        .*sqrt(factorial(m)*factorial(n)./factorial(m-k)./factorial(n-k))...
        .*O_mk_g.*O_g_nk;
end
%
% FIRST INDEX, MM, IS EXCITED STATE
% SECOND INDEX, NN, IS GROUND STATE
%
ovlp(mm,nn)=O_0_0^-1 * (num);
%
end
end
