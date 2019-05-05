%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matout] = ndaug(matin,n)

% Stacks along the main diagonal "n" copies of "matin".
%
% USAGE  : [matout] = ndaug(matin,n)
% INPUT  : matin    = a SYSTEM/VARYING/CONSTANT matrix
%          n        = how many copies of "matin" will be diagonally augmented
% OUTPUT : matout   = resulting output SYSTEM/VARYING/CONSTANT matrix
%
% Alexander Lanzon - 5 May 2000

if (nargin == 0) | (nargin > 2)
  disp('USAGE: [matout] = ndaug(matin,n)');
  return;
end

accum = matin;
for kk = 2:n
  accum = daug(accum,matin);
end
matout = accum;

return;