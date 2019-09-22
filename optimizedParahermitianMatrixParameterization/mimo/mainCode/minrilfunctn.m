%
%Function to obtain the minimal realisation of state-space realization
%given the state space matrices 
%

function [a,b,c,d] = minrilfunctn(aminril,bminril,cminril,dminril)
testminrildummy = ss(aminril,bminril,cminril,dminril);
dummy = minreal(testminrildummy);
[a,b,c,d] = ssdata(dummy);
return;