function M=rcbasis(nquanta,nmodes,target)
% Recursively generate integer partition with sum constraint
%   nquanta         : sum constraint
%   nmodes          : number of integers in sum
%   target             : remaining number of quanta to assign

% Test for simple case when the number of columns to assign target number
% of quanta is only 1.
if nmodes==1
    % Check that there is not an error in how many total quanta need to be
    % assigned versus how many currently need to be assigned. target should
    % not be greater than nquanta
    if target<=nquanta
        M=target;
    else
        M=[];
    end
% Recursive algorithm for when there are multiple columns to assign nquanta
% 
else
    M=[];
    for ii=0:min(nquanta,target)
        v=rcbasis(nquanta,nmodes-1,target-ii);
        M=[M;[repmat(ii,size(v,1),1) v]];
    end
end
end
