function [fcall] = fcinfo_I3(base,disp,nmode,nquanta)
% Calculate product of overlap integrals for all basis states
%   base            : matrix of quantum number for each basis state
%   disp             : unitless displacements for vibrational modes
%   nmode         : numbr of vibrational modes
%   nquanta       : number of total quanta to compute

%% Compute and store overlap integrals for all modes
[nstates,~]=size(base);

% Loop over vibrational modes
ovlpall_tc=zeros(nquanta+1,nquanta+1,nmode,'double');
for in=1:nmode
    d=disp(in);
    % Temporary define ground and excited state frequencies to be equal.
    wn=1;wm=1;
    % Calculate overlap integral for mode in
    ovlpall_tc(:,:,in)=fcfac2_TC(d,nquanta,wn,wm);
end % End loop over vibrational modes
% Reorder dimensions of object holding overlap integrals.
ovlpall_tc=permute(ovlpall_tc,[3 1 2]);
      
%% Compute and store products of overlap integrals for all basis states
% Outer loop over number of states

fcall=zeros(nstates,nstates,'double');
for iqO=1:nstates
    % Inner loop over number of states
    for iqI=1:nstates
        % Set product equal to 1
        num=1;
        % Loop over vibrational modes
        iv=1:nmode;
            % Get first quantum number
            j1=base(iqO,iv);
            % Get second quantum numbr
            j2=base(iqI,iv);
            % Multiply overlap integrals
            temp=ovlpall_tc(iv,j1+1,j2+1);
            for ii=1:nmode
                num=num*temp(ii,ii,ii);
            end
        % End loop over vibrational modes
        % Store product of overlap integrals
        fcall(iqO,iqI)=num;
    end % End inner loop over number of states
end % End outer loop over number of states

% fcall=zeros(nstates,nstates,'double');
% for iqO=1:nstates
%     % Inner loop over number of states
%     for iqI=1:nstates
%         % Set product equal to 1
%         num=1;
%         % Loop over vibrational modes
%         for iv=1:nmode
%             % Get first quantum number
%             j1=base(iqO,iv);
%             % Get second quantum numbr
%             j2=base(iqI,iv);
%             % Multiply overlap integrals
%             num=num*ovlpall_tc(iv,j1+1,j2+1);
%         end % End loop over vibrational modes
%         % Store product of overlap integrals
%         fcall(iqO,iqI)=num;
%     end % End inner loop over number of states
% end % End outer loop over number of states
