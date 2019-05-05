    %
    %Program that would take a system matrix and approximate the
    %state-space matrices to a particular user determined decimal places
    %
   
    function [matout] = roundupsys(matin,k)
    %
    %Decompose the input matrix into its transfer function form
    %
    [Adummy,Bdummy,Cdummy,Ddummy] = unpck(matin);
                            %
                            %Approximate to a predetermined decimal places the state-space matrix of M
                            %
                            %
                            Adummydp = Adummy;
                            [rnum,cnum] = size(Adummydp);
                                rcnt = 1;
                                while rcnt <= rnum
                                    ccnt = 1;
                                    while ccnt <= cnum
                                    Adummydp(rcnt,ccnt) = (round(Adummydp(rcnt,ccnt)*10^k))/10^k;
                                    ccnt = ccnt + 1;
                                    end
                                    rcnt = rcnt + 1;
                                end


                            %
                            Bdummydp = Bdummy;
                            [rnum,cnum] = size(Bdummydp);
                                rcnt = 1;
                                while rcnt <= rnum
                                    ccnt = 1;
                                    while ccnt <= cnum
                                    Bdummydp(rcnt,ccnt) = (round(Bdummydp(rcnt,ccnt)*10^k))/10^k;
                                    ccnt = ccnt + 1;
                                    end
                                    rcnt = rcnt + 1;
                                end

                            %
                            Cdummydp = Cdummy;
                            [rnum,cnum] = size(Cdummydp);
                                rcnt = 1;
                                while rcnt <= rnum
                                    ccnt = 1;
                                    while ccnt <= cnum
                                    Cdummydp(rcnt,ccnt) = (round(Cdummydp(rcnt,ccnt)*10^k))/10^k;
                                    ccnt = ccnt + 1;
                                    end
                                    rcnt = rcnt + 1;
                                end
                            %
                            Ddummydp = Ddummy;
                            [rnum,cnum] = size(Ddummydp);
                                rcnt = 1;
                                while rcnt <= rnum
                                    ccnt = 1;
                                    while ccnt <= cnum
                                    Ddummydp(rcnt,ccnt) = (round(Ddummydp(rcnt,ccnt)*10^k))/10^k;
                                    ccnt = ccnt + 1;
                                    end
                                    rcnt = rcnt + 1;
                                end
                                
   %repacking the system matrix and getting it ready for output 
   matout = pck(Adummydp,Bdummydp,Cdummydp,Ddummydp);
    
   return;