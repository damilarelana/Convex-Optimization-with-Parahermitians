    %
    %Program that would take a system matrix and multiplies the C and D matrix by a
    %positive scalar value k 
    %
   
    function [matout] = syscalarmult(matin,k)
    %
    %Decompose the input matrix into its transfer function form
    %
    [Adummy,Bdummy,Cdummy,Ddummy] = unpck(matin);
    Cdummy = Cdummy*(10^k);
    Ddummy = Ddummy*(10^k);
                                
   %repacking the system matrix and getting it ready for output 
   matout = pck(Adummy,Bdummy,Cdummy,Ddummy);
    
   return;