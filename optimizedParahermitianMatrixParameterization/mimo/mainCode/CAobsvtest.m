 function [flagctrdummy] = CAobsvtest(Adummy,Cdummy)
    %
    flagctrdummy = 'Blank';
   %
   %Obtain the eigenvalues
   %
   %
   %Start preparing the extracted positive real eigenvalues for testing via
   
   %
   %Obtain the controllability matrix
   %for controllability
   %
   Obmatrix = obsv(Adummy,Cdummy);
   
   [rdummyco,cdummyco] = size(Obmatrix);
      
   %
   rankobs = rank(Obmatrix);
       if (rankobs == cdummyco)
           flagctrdummy = 'Yes';
       else
           flagctrdummy = 'No';
       end
   
   return;