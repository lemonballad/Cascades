clear;
% close all;
%
%
%
gam_eg=2010;
for j=1:10
d=(j-1)*0.5+0.1;
displacement(j)=d;
[ratio] = cascadesub4(d,gam_eg);
s(j)=ratio;
end
figure;semilogy(displacement,s);
% contourf(dephasing,displacement,log10(s));colorbar;

% dlmwrite('.\ratio.dat',s);dlmwrite('.\displacement.dat',displacement);
%*3.8865e147.7729e14