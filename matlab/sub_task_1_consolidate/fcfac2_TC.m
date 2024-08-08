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
O_0_0=exp(-s/2);

%% 
for mm=1:nquanta+1
for nn=1:nquanta+1
m=mm-1; %actual quantum number for excited state
n=nn-1; %actual quantum number for ground state
%
num=0;
for k=0:min(mm,nn)-1;
O_mk_g=O_0_0 * (-1).^(m-k) .* s.^((m-k)/2)...
    ./sqrt(factorial(m-k));

O_g_nk=O_0_0 * s.^((n-k)/2)...
    ./sqrt(factorial(n-k));

%     num=num+(2*sqrt(R)/(1+R)).^k...
%         .*1./factorial(k)...
%         .*sqrt(factorial(m)*factorial(n)./factorial(m-k)./factorial(n-k))...
%         .*O_mk_g.*O_g_nk;
    num=num+1 ...
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
