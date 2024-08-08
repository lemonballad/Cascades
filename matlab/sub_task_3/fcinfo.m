function [fcall] = fcinfo(base,disp,nmode)
% c
% c     COMPUTE AND STORE OVERLAP INTEGRALS FOR ALL MODES
% c
[nstates,ij]=size(base);
% ovlpall=zeros(1:20,1:20,1:20);
%
%
%
      for i=1:nmode
      d=disp(i);
      [fc]=fcfac2(d);
      [ilim,ilim]=size(fc);
% c
      for ii=1:ilim
      for iii=1:ilim
      ovlpall(i,ii,iii)=fc(ii,iii);
      end
      end
% c      
      end
% c
% c     COMPUTE AND STORE PRODUCTS OF OVERLAP INTEGRALS FOR ALL BASIS STATES
% c
      for i=1:nstates
      for ii=1:nstates
% c
      num=1;
      for iii=1:nmode
      j1=base(i,iii);
      j2=base(ii,iii);
      num=num*ovlpall(iii,j1+1,j2+1);
      end
% c
      fcall(i,ii)=num;
% c      
      end
      end
