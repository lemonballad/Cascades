function [base,wviball] = basis(nmode,wvib);
%
%
%
%
%
% c
% c     SPECIAL CASE FOR 1 MODE
% c
      if nmode==1 
      iq=2;
      ic=1;
      for i=0:iq
      for ii=0:iq
% c
      if ii==i
      base(ic,1)=ii;
      ic=ic+1;
      end
      nstates=ic-1;
% c            
      end
      end
      end
% c
% c     SPECIAL CASE FOR 2 NMODES
% c
      if nmode==2 
      iq=7
      ic=1
      for i=0:iq
      for ii=0:iq
      for iii=0:iq
% c
      if ii+iii==i
      base(ic,1)=ii;
      base(ic,2)=iii;
      ic=ic+1;
      end
      nstates=ic-1;
% c            
      end
      end
      end
      end
% c
% c     SPECIAL CASE FOR 3 MODES
% c
      if nmode==3 
      iq=7
      ic=1
      for i=0:iq
      for ii=0:iq
      for iii=0:iq
      for iiii=0:iq
% c
      if ii+iii+iiii==i 
      base(ic,1)=ii;
      base(ic,2)=iii;
      base(ic,3)=iiii;
      ic=ic+1;
      end
      nstates=ic-1;
% c            
      end
end
end
end
      end
% c
% c     COMPUTE ENERGIES OF BASIS STATES HERE
% c
      for i=1:nstates
      num=0;
      for ii=1:nmode
      j=base(i,ii);
      num=num+(j+0.5)*wvib(ii);
      end
      wviball(i)=num;
      end 