 function [flagctr] = ABcontrtest(Adummy,Bdummy)
    %
    flagctr = 'Blank';
   %
   %Obtain the eigenvalues
   %
   %
   %Start preparing the extracted positive real eigenvalues for testing via
   
   %
   %Obtain the controllability matrix
   %for controllability
   %
   Comatrix = ctrb(Adummy,Bdummy);
   
   [rdummyco,cdummyco] = size(Comatrix);
     
   %
   rankctr = rank(Comatrix);
       if (rankctr == rdummyco)
           flagctr = 'Yes';
       else
           flagctr = 'No';
       end
   
   return;