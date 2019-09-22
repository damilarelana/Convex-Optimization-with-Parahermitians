%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT] = nconv(IN,n)

% This function convolves input polynomial row vector 'IN' with itself
% 'n' times and outputs the result in variable 'OUT'.
% 
% USAGE  : [OUT] = nconv(IN,n)
% INPUT  : IN    = an input polynomial row vector
%          n     = how many times this vector will be convolved with itself
% OUTPUT : OUT   = resulting output polynomial row vector
%
% Alexander Lanzon - 15 May 2000

if (nargin > 2)
  disp('USAGE: [OUT] = nconv(IN,n)');
  return;
end

accum = IN;
for kk = 2:n
  accum = conv(accum,IN);
end
OUT = accum;

return;