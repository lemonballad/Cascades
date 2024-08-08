clear all
clc

ratio=0.25;% ratio=b_thickness/a_thickness
a=10;
b=a*(1+ratio);
n=1:50;
x=-b:0.1:b;
v0=-0.005;
E=0;
V=v0*ones(1,length(x));V(abs(x)>a)=0;
hbar=1;
m=1;
L=2*b;
for in=n
    pib(:,in)=sqrt(2/L)*cos(in*pi*x/L);
    d2pib(:,in)=-(L/in/pi)^2*sqrt(2/L)*cos(in*pi*x/L);
end
plot(x,V,x,pib(:,1),x,d2pib(:,1));

for in=n
    for ij=n
        Knj(:,in,ij)=(-hbar^2/2/m)*pib(:,in)*d2pib(ij);
        Vnj(:,in,ij)=pib(:,in)*V*pib(:,ij);
    end
end

H=Knj+Vnj;
for ix=1:length(x)
[c,e]=eig(permute(H(ix,:,:),[2 3 1]));
Cnj(ix,:,:)=c;Enj(ix,:)=abs(diag(e));
end
in=10;
plot(x,V,x,conj(Cnj(:,in,in)).*Cnj(:,in,in)/in+Enj(end,in));

