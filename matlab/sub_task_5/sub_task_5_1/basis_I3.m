function [base,wviball] = basis_I3(nmode,nquanta,wvib)
% Create basis set and energies
%   nmode       : number of vibrational modes included in basis set
%   nquanta     : cutoff for number of quanta in basis set
%   wvib           : set of vibrational modes in wavenumbers with length nmode
%
%   base           : returns a nstates by nmode matrix of integers
%   (nstates      : binomial coefficient ==> (nquanta+nmode)!/nquanta!/nmode!)
%   wviball        : array of vibrational energies for each basis state
%
% Example:
%  basis(2,2,w)
%   base = [0 0; 0 1; 1 0; 0 2; 1 1; 2 0]
%   wviball = 

%% Throw Errors if any
% Is nmode greater than number of modes in wvib
if nmode>length(wvib)
    error('Error:\n%s%s\n%s%d%s%d',...
        'The variable "nmode" corresponds to more vibrational modes'... 
        ,' than reported in "wvib".'...
        ,'nmode = ',nmode,' and wvib = ',length(wvib));
end

% Is wvib empty
if isempty(wvib)
    error('Error: The variable wvib is empty');
end

%% Generate set of basis states
base=[];
% Loop over number of quanta
for nq=0:nquanta
    base=[base;rcbasis(nq,nmode,nq)]; %#ok<AGROW>
end % End loop over number of quanta
% Calculate number of states in basis
nstates=nchoosek(nquanta+nmode,nmode);

%% Compute Energies for basis states
wviball=zeros(1,nstates,'double');
% Loop over basis states
for is=1:nstates
    % define vector of quanta for basis state
    j=base(is,:);
    % Calculate vibrational energy for basis state
    wviball(is)=sum((j+0.5).*wvib(1:nmode));
end % End loop over basis states