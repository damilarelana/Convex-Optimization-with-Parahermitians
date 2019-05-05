    %
    %Program that would take a system matrix and check if the realization
    %is detectable and stabilizable or not using the PBH test. It outputs the flag "Yes" indicating that its detectable and stabilizable 
    %
   
    function [flagstab,flagdect] = stabdecttest(matin)
    %
    flagstab = 'Blank';
    %
    flagdect = 'Blank';
    
    %
    %Decompose the input matrix into its transfer function form
    %
    [Adummy,Bdummy,Cdummy,Ddummy] = unpck(matin);
    [radummy,cadummy]= size(Adummy);
   
    %
   %Obtain the eigenvalues
   %
   Adummyeig = eig(Adummy);
   
   %
   %Extract the real parts of the eigen values
   %
   ReAdummyeig = real(Adummyeig);
   
   %
   %Determine size of the real part of the eigen value vector 
   %
   [rowdummy,coldummy] = size(ReAdummyeig);
   
   %Start the process of extracting the vectors that have positive real
   %parts
   pbheig = [];
   cnt = 1;
   while cnt <= rowdummy
       if (ReAdummyeig(cnt) >= 0) 
          pbheig = abv(pbheig,Adummyeig(cnt));
       end
       cnt = cnt + 1;
   end
   
   if (isempty(pbheig) == 0)
               %
               %Start preparing the extracted positive real eigenvalues for testing via
               %
               %for stabilizability
               [rdummypbh,cdummypbh] = size(pbheig);
               pbhcnt = 1;
               while pbhcnt <= rdummypbh
                   Atestdummy = Adummy - (pbheig(pbhcnt)*eye(coldummy));
                   Teststabmatrix = [Atestdummy Bdummy];
                   [rteststab,cteststab] = size(Teststabmatrix);
                   rankstab = rank(Teststabmatrix);
                   if (rankstab == rteststab)
                       flagstab = 'Yes';
                   else
                       flagstab = 'No';
                       break;
                   end
                   pbhcnt = pbhcnt + 1; 
               end
                %for detectability
               [rdummypbh,cdummypbh] = size(pbheig);
               pbhcnt = 1;
               while pbhcnt <= rdummypbh
                   Atestdummy = Adummy - (pbheig(pbhcnt)*eye(coldummy));
                   Testdectmatrix = abv(Atestdummy,Cdummy);
                   [rtestdect,ctestdect] = size(Testdectmatrix);
                   rankdect = rank(Testdectmatrix);
                   if (rankdect == ctestdect)
                       flagdect = 'Yes';
                   else
                       flagdect = 'No';
                       break;
                   end
                   pbhcnt = pbhcnt + 1; 
               end
   end   
   
   return;