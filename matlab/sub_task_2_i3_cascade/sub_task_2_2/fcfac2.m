function [ovlp] = fcfac2(d);
%
% COMPUTE OVERLAP INTEGRALS BY HELLER'S METHOD
% EQUATION 7 IN MYERS ET. AL. J. CHEM. PHYS., VOL=77, PAGE=3857, YEAR=1982.
%
s=.5*d^2;
R=2;
%
%
lim=10;
for mm=1:lim
for nn=1:lim
m=mm-1; %actual quantum number for excited state
n=nn-1; %actual quantum number for ground state
%
num=0;
klimit=mm;
if nn<mm
klimit=nn;
end
for kk=1:klimit
k=kk-1;
Omkg=exp(-s/2)*(-1)^(m-k)*s^((m-k)/2)/sqrt(factorial(m-k));
Ognk=exp(-s/2)*s^((n-k)/2)/sqrt(factorial(n-k));
num=num+(2*sqrt(R)/(1+R))^k*1/factorial(k)*sqrt(factorial(m)*factorial(n)/factorial(m-k)/factorial(n-k))...
    *Omkg*Ognk;
end
%
% FIRST INDEX, MM, IS EXCITED STATE
% SECOND INDEX, NN, IS GROUND STATE
%
ovlp(mm,nn)=exp(s/2)*num;
Omkg
Ognk
num
%
end
end
