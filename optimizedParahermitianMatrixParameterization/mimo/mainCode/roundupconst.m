    %
    %Program that would take a constant matrix and approximate the
    %to a particular user determined decimal places
    %
   
    function [matout] = roundupconst(matin,k)
    %Approximate to a predetermined decimal places the constant matrix
    %given
                            %
                            %
                            Constdummydp = matin;
                            [rnum,cnum] = size(Constdummydp);
                                rcnt = 1;
                                while rcnt <= rnum
                                    ccnt = 1;
                                    while ccnt <= cnum
                                    Constdummydp(rcnt,ccnt) = (round(Constdummydp(rcnt,ccnt)*10^k))/10^k;
                                    ccnt = ccnt + 1;
                                    end
                                    rcnt = rcnt + 1;
                                end


                            
                                
   %repacking the system matrix and getting it ready for output 
   matout = Constdummydp;
    
   return;