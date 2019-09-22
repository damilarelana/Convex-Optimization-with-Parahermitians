%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       CODE WRITTEN BY LANA DAMILARE YOMI (Copyright November 2009)     % 
%                           FOR HIS MPHIL THESIS 2009                    % 
%                               AT THE                                   %
%                         UNIVERSITY OF MANCHESTER                       %
%         (Note that this code also contains some useful functions       %
%              written by my Supervisor Dr. Alexander Lanzon)            %  
%                                                                        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete compareMsys.jpg
delete compareMsys.eps
delete wrkspcbfmu
delete compareds.jpg
delete compareds.eps
delete comparegammas.jpg
delete comparegammas.eps
%
%
close all
clear all
clear classes
clear java
clear wrkspcbfmu
clear myvariables
clear Dfull
clear Decon
clc
%
%
%Diary begins to take note of what is going on
%
diary wrkspcbfmu
disp(' --------------------------------------------------------- ')
disp(' ------------        DIARY BEGINS       ------------------ ')
disp(' --------------------------------------------------------- ')

%
%
%Publishing the results of the code execution
%publish('Latest','html')
%open html/Latest.html
%Getting necessary inputs
%
%
%Specifying the frequency response for the inputs
%
omega = logspace(-4,4,100);

%
tminlimit = -1;
%
%inputing the plant to be used
%
%
Aginpt = [-19.48 -0.91 -3.15 7.94 -4.57 -4.64; -1.02 -10.59 -0.64 23.82 -13.71 -13.92; -3.24 -0.70 -30.99 15.42 -8.88 -9.01; 0 -0.15 -0.12 -5.29 -6.24 -7.00; 0 2.14 1.72 14.97 -24.21 -14.59; 0 2.70 2.17 10.21 -11.55 -39.76];
Bginpt = [1.25 -0.38 -1.27 0 -0.02 0.76; -0.64 -0.19 0.98 0 0.62 0; 0.58 -1.38 -0.04 0.86 0.47 0; 0 0.03 0 0 -0.05 0; 0 -0.36 0 0 0.71 0; 0 -0.45 0 0 0.89 0];
Cginpt = [-23.77 8.66 4.30 13.44 -7.73 -7.85; -2.74 0.24 0.19 17.37 -10.00 -10.15; 0 -9.23 -24.44 0 0 0; 0 -0.54 -10.85 11.91 -6.85 -6.96; -5.11 1.44 11.20 -2.70 1.55 1.58; -0.02 0.16 5.95 11.73 -6.75 -6.86];
Dginpt = [4.14 8.55 0 9.24 0.61 -5.50; 0 -0.40 3.52 0 0.79 0.40; 0 -7.26 11.33 0 0 0; 3.18 -4.72 0 0 0.54 0; 0 -6.07 0 10.28 -0.12 0; 0 -0.27 -0.52 3.95 0.53 11.53];
Msys = pck(Aginpt,Bginpt,Cginpt,Dginpt);
                
%Obtain frequency response

Mfrq = frsp(Msys,omega); 


disp(' ')
disp(' ------------------------------------------------------------------ ')
disp(' ')

%
%Testing the controllability and observability of Msys(s)  
%
[flagctrm] = ABcontrtest(Aginpt,Bginpt);
if (strcmp(flagctrm,'Yes')== 1)
disp(' ')
disp('Msys(s) is controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
else 
disp(' ')
disp('Msys(s) is not controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
end

[flagobsrm] = CAobsvtest(Aginpt,Cginpt);
if (strcmp(flagobsrm,'Yes')== 1)
disp(' ')
disp('Msys(s) is observable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
else 
disp(' ')
disp('Msys(s) is not observable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
end

%
%For Msys
%
figure
[uMsys,sMsys,vMsys] = vsvd(Mfrq);
vplot('liv,lm',sMsys);
ylabel('Singular Value of  Msys(s)','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -depsc compareMsys
print -djpeg compareMsys
%
%
%

disp('The M lower lft matrix has ')
minfo(Msys)
disp('and its packed form is')
seesys(Msys,'%13.7f') 
disp(' ')
disp('and its eigenvalues are ')
eig(Aginpt)
disp(' ')
%
%
Mdpsys = roundupsys(Msys,2);
%
%
disp('The Mdpsys lower lft matrix is the approx of Msys and it has ')
minfo(Mdpsys)
disp('and its packed form is')
seesys(Mdpsys,'%13.2f') 
disp(' ')
%specifically done to crosscheck the approximation done using system matrix
[aMdpsys,bMdpsys,cMdpsys,dMdpsys] = unpck(Msys);
aMdpsys = roundupconst(aMdpsys,2);
disp('aMdpsys ')
seesys(aMdpsys,'%13.2f')
bMdpsys = roundupconst(bMdpsys,2);
disp('bMdpsys ')
seesys(bMdpsys,'%13.2f')
cMdpsys = roundupconst(cMdpsys,2);
disp('cMdpsys ')
seesys(cMdpsys,'%13.2f')
dMdpsys = roundupconst(dMdpsys,2);
disp('dMdpsys ')
seesys(dMdpsys,'%13.2f')
disp(' --------------------------------------------------------- ')


%
%Clear the variables to free up memory space
%
clear vMsys
clear sMsys
clear uMsys
clear aMdpsys
clear bMdpsys
clear cMdpsys
clear dMdpsys

%
%
%Obtaining the dimension of the Lower LFT, Please note that from now
%onwards that Msys represents the choice of the user i.e. it can either
%have been a suboptimal or an optimal controller
%
[Mtype,Mrow,Mcol,Mnum] = minfo(Msys);

%
%Extract the state-space matrices of
%Msys to be used in the LMI program
%
[Ag,Bg,Cg,Dg] = unpck(Msys);

%
%Extract the dimensions of the state matrices of Msys
%

[agrow,agcol] = size(Ag);
[bgrow,bgcol] = size(Bg);
[cgrow,cgcol] = size(Cg);
[dgrow,dgcol] = size(Dg);

%
%Initialise the while loop breakers to be used to break from the tau and N
%while loops once there is dimension incompatibility between Laguerre state
%matrices and the original plant M(s)
%
nonrecurmchktauNbrk = 0;
%
%Obtain the uncertainty structure ready for usage
%
deltatoscalar = 2;
scalartodscales = 4;

%-----------------------------------------------------------------
               
                              
                %
                %
                % Non-recursive method starts here
                %
                %
                        Nchkcnt = 0;
                        while Nchkcnt < 1
                        disp(' ')
                        disp(' ')
                        N = input('What is your positive integer value of N to be used in the Laguerre Parameterisation? ');
                        disp(' ')
                                if (N >= 1)
                                   Nchkcnt = 1;
                                else
                                disp(' ')
                                disp('Please make sure N is an positive integer value greater than 1')
                                disp(' ')
                                disp(' ')  
                                end
                        end
                        %
                        tauchkcnt = 0;
                        while tauchkcnt < 1
                        disp(' ')
                        disp(' ')
                        tau = input('What is your positive and real value of \tau to be used in the Laguerre Parameterisation? ');
                                if ((isreal(tau) == 1) && (tau > 0))
                                tauchkcnt = 1;
                                else
                                disp(' ')
                                disp(' Only positive and real values of tau is allowed')
                                disp(' ')
                                disp(' ')
                                end
                        end        



                                    %--------------------------------------------------------------------------
                                    %
                                    %Synthesis the initial A_d and B_d from the B(s) of the Laguerre Parameterisation using the latest version made today 9th September 2009
                                    %Please note that for now this program can only handle cases for which the Dscale has only 1 symetric scalar part and 1 symetric full matrix part on its
                                    %diagonal.
                                    %--------------------------------------------------------------------------
                                    %Start constructing Ad and Bd from the
                                    %corresponding Ad11,Ad22, B11 and B22
                                    %using q1 = deltatoscalar and q2 = scalartodscales
                                    %
                                    %Assign values with respect to the
                                    %defined uncertainty structure
                                    %
                                    q1 = deltatoscalar; 
                                    q2 = scalartodscales;
                                    
                                    %
                                    %
                                    %Part to obtain the scalar part of the
                                    %Dscale structure note that its
                                    %transpose just like laguerre
                                    %Note that since the dscale is supposed
                                    %Here the simple SISO plant d1(s) is
                                    %approximated in a way that gives a 1 by 1 dimension and then distributed
                                    %across the diagonal q1 times i.e.
                                    %since d1(s) occurs q1 times on the
                                    %diagonal.
                                    %
                                  
                                    %
                                    %Laguerre Implementation
                                    %
                                                                        
                                    numdummy = 1;
                                    dendummy = nconv([1 2/(tau)],N);
                                    Hdummy = tf(numdummy,dendummy);
                                                                        
                                    %
                                    [a1dummy,b1dummy,c1dummy,d1dummy] = tf2ss(1,nconv([1 2/(tau)],N));
                                    
                                    %
                                    %Before spreading along the diagonal
                                    %
                                    
                                    disp(' ')
                                    disp('Before spreading using Iq, we have that Laguerre is ')
                                    disp(' ')
                                    minfo(a1dummy)
                                    seesys(a1dummy,'%13.2f') 
                                    disp(' ')
                                    disp('and its eigenvalues are ')
                                    eig(a1dummy)
                                    disp(' ')
                                    
                                    %
                                    %Check whether or not the A and B
                                    %matrix are controllable or not after
                                    %ndauging
                                    %
                                    [flagctra] = ABcontrtest(a1dummy,b1dummy);
                                    
                                    if (strcmp(flagctra,'Yes')== 1)
                                    disp(' ')
                                    disp('(sisoA,sisoB) is controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    elseif (strcmp(flagctra,'No')== 1)
                                    disp(' ')
                                    disp('(sisoA,sisoB) is not controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    else
                                    disp(' ')
                                    disp('Controllability of (sisoA,sisoB) is unknown')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    end
                                   
                                    %
                                    %
                                    %
                        %
                                    %FOR AD1
                                    %                                   
                                                                                                
                                    %Spreading along the diagonal in laguerre
                                    
                                    a1dummy2B1mnrl = ndaug(a1dummy,q1);
                                    b1dummy2B1mnrl = ndaug(b1dummy,q1);
                                                                        
                                    disp(' ')
                                    disp('After spreading using Iq1, we have that MIMO Ad1 is ')
                                    disp(' ')
                                    minfo(a1dummy2B1mnrl)
                                    seesys(a1dummy2B1mnrl,'%13.2f') 
                                    disp(' ')
                                    disp('and its eigenvalues are ')
                                    eig(a1dummy2B1mnrl)
                                    disp(' ')
                                    
                                    %
                                    %Check whether or not the A and B
                                    %matrix are controllable or not after
                                    %ndauging
                                    %
                                    [flagctrb] = ABcontrtest(a1dummy2B1mnrl,b1dummy2B1mnrl);
                                     
                                    if (strcmp(flagctrb,'Yes')== 1)
                                    disp(' ')
                                    disp('(Ad1,Bd1) is controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    elseif (strcmp(flagctrb,'No')== 1)
                                    disp(' ')
                                    disp('(Ad1,Bd1) is not controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    else
                                    disp(' ')
                                    disp('Controllability of (Ad1,Bd1) is unknown')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    end
                                    
                                    
                                    %--------------------------------------
                                    %Note that in contrast to that of
                                    %Bs(s), we need to construct the whole
                                    %B1(s) first with the one before
                                    %scaling
                                    %Now attempting to form the B1(s) as obtained in the laguerre parameterization                                    %
                                    %--------------------------------------
                                    %Augment Dscales so as to be able to produce Iq1 within B1(s)
                                    [rowsclra,colsclra]= size(a1dummy2B1mnrl);
                                    [rowsclrb,colsclrb]= size(b1dummy2B1mnrl);

                                    %
                                    %Scale to form the big scalar diagonal version
                                    %
                                    Ascalarbig = a1dummy2B1mnrl;
                                    [rowAscalarbig,colAscalarbig] = size(Ascalarbig);
                                    %
                                    %assign value to use in the LMI
                                    %
                                    SizeofAd1 = rowAscalarbig;

                                    %
                                    Bscalarbig = b1dummy2B1mnrl;
                                    [rowBscalarbig,colBscalarbig] = size(Bscalarbig);
                                    %
                                    
                                    %
                                    %FOR AD2
                                    %
                                    
                                    %Spreading along the diagonal in laguerre
                                 
                                    a2dummy2B2mnrl = ndaug(a1dummy,q2);
                                    b2dummy2B2mnrl = ndaug(b1dummy,q2);
                                                                        
                                    disp(' ')
                                    disp('After spreading using Iq2, we have that MIMO Ad2 is ')
                                    disp(' ')
                                    minfo(a2dummy2B2mnrl)
                                    seesys(a2dummy2B2mnrl,'%13.2f') 
                                    disp(' ')
                                    disp('and its eigenvalues are ')
                                    eig(a2dummy2B2mnrl)
                                    disp(' ')
                                    
                                    %
                                    %Check whether or not the A and B
                                    %matrix are controllable or not after
                                    %ndauging
                                    %
                                    [flagctrd] = ABcontrtest(a2dummy2B2mnrl,b2dummy2B2mnrl);
                                    
                                    if (strcmp(flagctrd,'Yes')== 1)
                                    disp(' ')
                                    disp('(Ad2,Bd2) is controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    elseif (strcmp(flagctrd,'No')== 1)
                                    disp(' ')
                                    disp('(Ad2,Bd2) is not controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    else
                                    disp(' ')
                                    disp('Controllability of (Ad2,Bd2) is unknown')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    end
                                 
                                    %
                                    %
                                    [rowa2dummy2B2mnrl,cola2dummy2B2mnrl] = size(a2dummy2B2mnrl);
                                    %
                                    %assign value to use in the LMI
                                    %
                                    SizeofAd2 = rowa2dummy2B2mnrl;
                                    %
                                    %
                                    %
                                    [rowb2dummy2B2mnrl,colb2dummy2B2mnrl] = size(b2dummy2B2mnrl);


                                    

                                    %--------------------------------------
                                    %
                                    %CONSTRUCTING THE FINAL Ad and Bd THAT
                                    %WOULD BE USED WITHIN THE OPTIMISATION PROBLEM, 
                                    %such that A1 and A2 have the desired structure relative to the defined uncertainty structure
                                    %
                                    %--------------------------------------
                                   %
                                    %FOR Ad
                                    %
                                    
                                    Ad = daug(a1dummy2B1mnrl,a2dummy2B2mnrl);
                                    
                                    %
                                    disp(' ')
                                    disp('Ad is ')
                                    disp(' ')
                                    minfo(Ad)
                                    seesys(Ad,'%13.2f') 
                                    disp(' ')
                                    disp('and its eigenvalues are ')
                                    eig(Ad)
                                    disp(' ')
                                    
                                    %
                                    %FOR Bd
                                    %
                                    
                                    Bd = daug(b1dummy2B1mnrl,b2dummy2B2mnrl);
                                    
                                    %
                                    disp(' ')
                                    disp('Bd is ')
                                    disp(' ')
                                    minfo(Bd)
                                    seesys(Bd,'%13.2f') 
                                    disp(' ')
                                    
                                    %
                                    %Check whether or not the A and B
                                    %matrix are controllable or not after
                                    %ndauging
                                    %
                                    [flagctrab] = ABcontrtest(Ad,Bd);
                                    
                                    if (strcmp(flagctrab,'Yes')== 1)
                                    disp(' ')
                                    disp('(Ad,Bd) is controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    elseif (strcmp(flagctrab,'No')== 1)
                                    disp(' ')
                                    disp('(Ad,Bd) is not controllable')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    else
                                    disp(' ')
                                    disp('Controllability of (Ad,Bd) is unknown')
                                    disp(' ')
                                    disp(' +++++++++++++++++++++++++++++ ')
                                    end
                                    
                                    %
                                    %Obtain the final large Ad and Bd
                                    %dimensions to be used in the LMI
                                    %programs
                                    %
                                    [Adrow,Adcol] = size(Ad);
                                    [Bdrow,Bdcol] = size(Bd);

                                    
                                    %
                                    %--------------------------------------------------------------
                                    %
                                    %
                                    %Factor to be used in the determination of the structure of
                                    %P in relation to the structure of Delta_hat
                                    %
                                    %For matrix to scalar using q1
                                    %
                                    [tvra1,tvca1]=size(Ascalarbig);
                                    if rem(tvra1,q1)== 0
                                    dummyfactor = tvra1/q1;
                                    disp(' ')
                                    disp([' The row of A-matrix for q1 is ',int2str(dummyfactor),' times the multiple of q1'])
                                    %
                                    sclrfactorfromq1 = dummyfactor;
                                    else
                                    disp(' ')
                                    disp(' Error: The row of  A-matrix for q1 is not an integer multiple of q1')
                                    %
                                    sclrfactorfromq1 = 0;
                                    end
                                    %
                                    %For scalar to matrix using q2
                                    %
                                    [tvra2,tvca2]=size(a2dummy2B2mnrl);
                                    if rem(tvra2,q2)== 0
                                    dummyfactor = tvra2/q2;
                                    disp(' ')
                                    disp([' The row of  A-matrix for q2 is ',int2str(dummyfactor),' times the multiple of q2'])
                                    %
                                    sclrfactorfromq2 = dummyfactor;
                                    else
                                    disp(' ')
                                    disp(' Error: The row of  A-matrix for q2 is not an integer multiple of q2')
                                    %
                                    sclrfactorfromq2 = 0;
                                    end

                                    %
                                    %------------------------------------------------------------------------
                                    %
                                    %


                                    %
                                    %------------------------------------------------------------------------
                                    %
                                    %

                                    %
                                    %Ensure that there is still row and
                                    %column compatibility between M(s) and
                                    %overall Bd and DdBs
                                    %

                                    %minimum version
                                    [btemprow,btempcol] = size(Bd);
                                    
                                    
                                    %
                                    %Comparison stage 
                                    %
                                    if (Mrow == btempcol)


                                    % Note that because i am dealing with a MIMO system, then the controller
                                    % canonical form is not advisable as it only complicates things further
                                    % Extracting the A and B matrices from B(s) dummy and constructing another
                                    % B(s) matrix that [I;0] output C matrix and [0;I] output D matrix
                                    %
                                    %
                                    [Adrow,Adcol] = size(Ad);
                                    [Bdrow,Bdcol] = size(Bd);
                                    %
                                    %Obtaining the approximated form 
                                    %
                                    Adapprox = roundupconst(Ad,3);
                                    Bdapprox = roundupconst(Bd,3);
                                    
                                    %
                                    %Form the big matrices Aehat and Behat to be used
                                    %during the LMI program
                                    %
                                  
                                    Aehat = [Ad Bd*Cg zeros(Adrow);zeros(agrow,Adrow) Ag zeros(agrow,Adrow);zeros(Adrow) zeros(Adrow,agrow) Ad];
                                    Behat = [Bd*Dg ; Bg ; Bd];
                                    
                                    [rplus2n,rplus2n] = size(Aehat);
                                    
                                    
                                                                        
                                    %
                                    %counter to ensure that programming module starts or does not start
                                    %
                                    nextcellcnt = 1;
                                    else
                                        disp(' ')
                                        disp(' ')
                                        disp('There is dimension imcompatibility between the row size of M(s) and the column size of the minimal realised Bd and Dd Laguerre state matrices')
                                        disp(' ')
                                        disp(' PROGRAM TERMINATED')
                                        nonrecurmchktauNbrk = 1;
                                        break
                                    end
                                    
                                    %
                                    %Clear the variables to free up memory space
                                    %
                                    clear dummyfactor
                                    clear Hdummy
                                    clear a1dummy
                                    clear b1dummy
                                    clear c1dummy
                                    clear d1dummy
                                  
                                    
                                    %%
                                    %
                                    %
                                    %--------- Start of Algorithm for Full Parameterization
                                    %
                                    %
                                    %obtain the tolerance for gammachoice
                                    %
                                    %
                                    gammaminref = eps;
                                    bsctnchkcnt = 0;
                                                    while bsctnchkcnt < 1
                                                        disp(' ')
                                                        disp(' ')
                                                        gammainputref = input('What is your choice of initial value of gamma? ');
                                                        disp(' ')
                                                        if (gammainputref > abs(eps)) 
                                                           bsctnchkcnt = 1;
                                                        else
                                                        disp(' ')
                                                        disp('Please only positive real values are allowed')
                                                        disp(' ')
                                                        algoquestion = input('Do you want to start with a default value of gamma ? Y/N: ','s');
                                                                if (strncmpi('y', algoquestion,1) == 1) || (strncmpi('yes', algoquestion,3) == 1)
                                                                    disp(' ')
                                                                    disp([' You have choosen to use default gamma = ',num2str(gammaminref),])
                                                                    disp(' ')
                                                                    gammainputref = gammaminref;
                                                                    bsctnchkcnt = 1;
                                                                end
                                                        end
                                                    end
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %                                                                        %  
                                    %         Now starting to write the main LMI program                     %
                                    %                                                                        %
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %
                                    %unpacking into its constituent system matrices
                                    %
                                    [Ai,Bi,Ci,Di] = unpck(Msys);
                                    %
                                    %Obtaining the "r" dimension while nothing that Ai is rxr, and Airow = r 
                                    %
                                    [Airow,Aicol] = size(Ai);
                                    %
                                    %
                                    %--------------------------------------------------------------------------
                                    %
                                    % First LMI with Full Decision Variables in the parahermitian
                                    %
                                    %--------------------------------------------------------------------------
                                    %

                                    %

                                    gammachoice = gammainputref;
                                    fullbisectioncnt = 0;
                                    xfeasfull = [];
                                    spctralcntfull = 0;
                                    %
                                    % Timer starts here
                                    tic
                                    %
                                    % Incremental Algorithm starts here
                                    %
                                    
                                    

                                        while gammachoice <= 100
                                              if ((spctralcntfull == 1) && (isempty(xfeasfull) == 0));
                                                  break
                                              else
                                                  spctralcntfull = 0;
                                                  xfeasfull = [];
                                              end
                                        fullbisectioncnt = fullbisectioncnt + 1;   
                                        
                                        % Now for the real lmi declaration noting that Mrow is the same as Bdrow
                                        %
                                        %Declare the create the varying matrix that var2con would then be used to
                                        %select data from within the LMI solution
                                        %
                                        %
                                        %
                                       
                                        disp([int2str(fullbisectioncnt),' Iteration of the Full LMI optimization using Gammachoice ',num2str(gammachoice),', N = ', int2str(N),' and tau = ',num2str(tau)])
                                       
                                        %
                                        %Declare the constant to use in the
                                        %LMI
                                        %
                                        Roneone = [eye(Adrow) zeros(Adrow,agrow) zeros(Adrow); zeros(Mrow,Adrow) Cg zeros(Mrow,Adrow); zeros(Adrow) zeros(Adrow,agrow) gammachoice*eye(Adrow)];
                                        Ronetwo = [zeros(Adrow,Mrow); Dg; zeros(Adrow,Mrow)];
                                        Rtwoone = [zeros(Mrow,Adrow) zeros(Mrow,agrow) zeros(Mrow,Adrow)];
                                        Rtwotwo = gammachoice*eye(Mrow);
                                        Routrfctor = [Roneone Ronetwo;Rtwoone Rtwotwo];
                                        
                                                                                
                                        %
                                        setlmis([])

                                        %setting up the main variables
                                        %nothing that Adrow = n, Mrow = m,
                                        %agrow = r Also note that Z reps -Z
                                        [Pdummy11,nPdummy11,sPdummy11] = lmivar(1,[SizeofAd1 1]); % for commuting property due to the Delta uncertainty stage two

                                        [Pdummy22,nPdummy22,sPdummy22] = lmivar(1,[SizeofAd2 1]); % for commuting property due to the scalar uncertainty stage two

                                        [P,nP,sP] = lmivar(3,[sPdummy11 zeros(SizeofAd1,SizeofAd2); zeros(SizeofAd2,SizeofAd1) sPdummy22]); %full n x n symmetric matrix
                                          
                                        [E,nE,sE] = lmivar(2,[Adrow Mrow]); %rectangular matrix E represents S an n x m matrix but replaced with E to avoid alphabetic confusion
                                     
                                        [Etranspose,nEtranspose,sEtranspose] = lmivar(3,sE.'); %rectangular matrix E represents S an n x m matrix but replaced with E to avoid alphabetic confusion
                                        
                                        [Rdummy11,nRdummy11,sRdummy11] = lmivar(1,[deltatoscalar 1]);
                                        
                                        [Rdummy22,nRdummy22,sRdummy22] = lmivar(1,[scalartodscales 1]);
                                        
                                        [R,nR,sR] = lmivar(3,[sRdummy11 zeros(deltatoscalar,scalartodscales); zeros(scalartodscales,deltatoscalar) sRdummy22]); %full m x m symmetric matrix
                                                                            
                                        [Z,nZ,sZ] = lmivar(1,[rplus2n 1]); %full n x n symmetric matrix
                                                                             
                                        [J,nJ,sJ] = lmivar(1,[Adrow 1]); %full n x n symmetric matrix
                                                                             
                                        [Drvoneone,nDrvoneone,sDrvoneone] = lmivar(3,[sP sE zeros(Adrow); sEtranspose sR zeros(Mrow,Adrow); zeros(Adrow) zeros(Adrow,Mrow) -1*sP]);

                                        [Drvonetwo,nDrvonetwo,sDrvonetwo] = lmivar(3,[zeros(Adrow,Mrow); zeros(Mrow); -1*sE]);
                                        
                                        [Drvtwoone,nDrvtwoone,sDrvtwoone] = lmivar(3,[zeros(Mrow,Adrow) zeros(Mrow) -1*sEtranspose]);
                                        
                                        [Drvtwotwo,nDrvtwotwo,sDrvtwotwo] = lmivar(3,-1*sR);
                                        %
                                        %likewise R has no imaginary component
                                        %
                                        %
                                        PRA = newlmi;

                                        %Declaring the lmi terms
                                

                                        %Structured singular value issue starts here
                                        lmiterm([1 0 0 0],Routrfctor) %
                                        lmiterm([1 1 1 Drvoneone],1,1) %
                                        lmiterm([-1 1 1 Z],Aehat.',1,'s')
                                        
                                        lmiterm([1 1 2 Drvonetwo],1,1) %
                                        lmiterm([-1 1 2 Z],1,Behat) %
                                        
                                        lmiterm([1 2 2 Drvtwotwo],1,1) %                                        

                                        %Positive definitedness of
                                        %parahermitian issue

                                        lmiterm([2 1 1 J],Ad.',1,'s') %
                                        lmiterm([-2 1 1 P],1,1) %
                                        lmiterm([-2 1 2 E],1,1) %
                                        lmiterm([2 1 2 J],1,Bd) %
                                        lmiterm([-2 2 2 R],1,1) %

                                        %positive definiteness
                                        %condition on the matrix Rspr
                                        %
                                        lmiterm([-3 1 1 R],1,1) %for full Rspr > 0
                                        
                                        

                                        
                                        %
                                        lmisysfull = getlmis;

                                        %Determining the LMI feasibility solution and assigning it to a matrix
                                        %storage
                                        [tminfull,xfeasfull] = feasp(lmisysfull,[0,200,0,25,0],tminlimit);
                                        if ((isempty(xfeasfull) == 0) && (tminfull < 0))
                                        % Piecewise extraction of the decision variables
                                        Poptfullrsltdummy = dec2mat(lmisysfull,xfeasfull,P);
                                        Soptfullrsltdummy = dec2mat(lmisysfull,xfeasfull,E);
                                        Roptfullrsltdummy = dec2mat(lmisysfull,xfeasfull,R);
                                        numfulldecvr = decnbr(lmisysfull);

                                        %
                                        LargeZdummy = dec2mat(lmisysfull,xfeasfull,Z);
                                        Jdummy = dec2mat(lmisysfull,xfeasfull,J);
%
                                        %--------------------------------------------------------------------------
                                        %Defining riccati equation parameters
                                        %(Ad-Bd*inv(R)*S')'*Y1 + Y1*(Ad-Bd*inv(R)*S') - Y1*Bd*inv(R)*Bd'*Y1 + P - S*inv(R)*S'
                                        % with P obtained as a full symetric matrix from the related optimisation
                                        % earlier done.
                                        %--------------------------------------------------------------------------
                                        %
                                                        %
                                                        Poptfullrslt = Poptfullrsltdummy;
                                                        Soptfullrslt = Soptfullrsltdummy;
                                                        Roptfullrslt = Roptfullrsltdummy;
                                                        %
                                                        %
                                                        Ahoptfullrslt = Ad - Bd*inv(Roptfullrslt)*Soptfullrslt';
                                                        Bhoptfullrslt = Bd*inv(Roptfullrslt)*Bd';
                                                        Qhoptfullrslt = Poptfullrslt - Soptfullrslt*inv(Roptfullrslt)*Soptfullrslt';

                                                        %Form the Hamiltonian matrix
                                                        Hammatrixfull = [Ahoptfullrslt -Bhoptfullrslt;-Qhoptfullrslt -Ahoptfullrslt'];
                                                        %
                                                        %Obtaining the stabilising solution
                                                        [Y11,Y12,failfull,reig_minfull] = ric_schr(Hammatrixfull);
                                                        if failfull == 0
                                                            Y1 = Y12*inv(Y11);
                                                            %
                                                        %
                                                        %increment counter
                                                        spctralcntfull = 1; %to show that this part of the loop was successful
                                                        %
                                                        else
                                                        disp('Problem with the determination of the riccati stabilising solution as per Full parameterisation')
                                                        end
                                        else
                                        disp('Problem with non-existence of a feasible solution as per Full parameterisation')
                                        end
                                        %
                                        %  
                                           Gammausednonrecursivefull = gammachoice;
                                           gammachoice = gammachoice + 0.1;
                                        %
                                        %
                                        end
                                        %
                                        %
                                        %
                                        if ((spctralcntfull == 1) && (isempty(xfeasfull) == 0))
                                        %
                                        %
                                        %Stop timer in the existence of a
                                        %solution 
                                        Fulltimer = toc;
                                        %
                                        %
                                        Nfull = N;
                                        taufull = tau;
                                        %
                                        %
                                        %
                                        disp(' ')
                                        disp(' --------------------------------------------------------------------------- ')
                                        disp(' ')
                                        disp(['Non-recursive Full opt Feasible with gamma = ',num2str( Gammausednonrecursivefull), ' , \tau = ', num2str(taufull), ' and N = ', int2str(Nfull)])
                                        disp(' ')
                                        disp(' --------------------------------------------------------------------------- ')
                                        disp(' ')
                                        %
                                        %
                                        %
                                        else
                                        %
                                        %Stop timer in case no solution was
                                        %obtained
                                        %
                                        Fulldummytimer = toc;
                                        %
                                        %
                                        %
                                        end
                                        
                                        
                                        %
                                        %--------------------------------------------------------------------------
                                        %
                                        % Second LMI with economised Decision Variables in the parahermitian
                                        %
                                        %--------------------------------------------------------------------------
                                        %
                                        %

                                        gammachoiceb = gammainputref;
                                        econbisectioncnt = 0;
                                        xfeasecon = [];
                                        spctralcntecon = 0;
                                        %
                                        % Timer starts here
                                        tic
                                        %
                                        % Algorithm starts here
                                        %
                                        while gammachoiceb <= 100
                                              if ((spctralcntecon == 1) && (isempty(xfeasecon) == 0));
                                                  break
                                              else
                                                  spctralcntecon = 0;
                                                  xfeasecon = [];
                                              end
                                        econbisectioncnt = econbisectioncnt + 1;   
                                        %
                                        %
                                        %
                                        
                                        disp([int2str(econbisectioncnt),' Iteration of the Economised LMI optimization using Gammachoice ',num2str(gammachoiceb),', N = ', int2str(N),' and tau = ',num2str(tau)])
                                        %
                                        %Declare the constant to use in the
                                        %LMI
                                        %
                                        Roneone = [eye(Adrow) zeros(Adrow,agrow) zeros(Adrow); zeros(Mrow,Adrow) Cg zeros(Mrow,Adrow); zeros(Adrow) zeros(Adrow,agrow) gammachoiceb*eye(Adrow)];
                                        Ronetwo = [zeros(Adrow,Mrow); Dg; zeros(Adrow,Mrow)];
                                        Rtwoone = [zeros(Mrow,Adrow) zeros(Mrow,agrow) zeros(Mrow,Adrow)];
                                        Rtwotwo = gammachoiceb*eye(Mrow);
                                        Routrfctor = [Roneone Ronetwo;Rtwoone Rtwotwo];
                                        
                                        %
                                        %
                                        %
                                        setlmis([])
                                        %setting up the main variables
                                        [Eecon,nEecon,sEecon] = lmivar(2,[Adrow Mrow]); %rectangular matrix E represents S an n x m matrix but replaced with E to avoid alphabetic confusion

                                        [Etransposeecon,nEtransposeecon,sEtransposeecon] = lmivar(3,sEecon.'); %rectangular matrix E represents S an n x m matrix but replaced with E to avoid alphabetic confusion

                                        [Recondummy11,nRecondummy11,sRecondummy11] = lmivar(1,[deltatoscalar 1]);

                                        [Recondummy22,nRecondummy22,sRecondummy22] = lmivar(1,[scalartodscales 1]);

                                        [Recon,nRecon,sRecon] = lmivar(3,[sRecondummy11 zeros(deltatoscalar,scalartodscales); zeros(scalartodscales,deltatoscalar) sRecondummy22]); %full m x m symmetric matrix

                                        [V,nV,sV] = lmivar(1,[Adrow 1]); %a full n x n symmetric matrix used for the showing that phi(s)'[0 S;S' R]phi(s) > 0 

                                        [Zecon,nZecon,sZecon] = lmivar(1,[rplus2n 1]); %full n x n symmetric matrix

                                        [Drvoneoneecon,nDrvoneoneecon,sDrvoneoneecon] = lmivar(3,[zeros(Adrow) sEecon zeros(Adrow); sEtransposeecon sRecon zeros(Mrow,Adrow); zeros(Adrow) zeros(Adrow,Mrow) zeros(Adrow)]);

                                        [Drvonetwoecon,nDrvonetwoecon,sDrvonetwoecon] = lmivar(3,[zeros(Adrow,Mrow); zeros(Mrow); -1*sEecon]);

                                        [Drvtwooneecon,nDrvtwooneecon,sDrvtwooneecon] = lmivar(3,[zeros(Mrow,Adrow) zeros(Mrow) -1*sEtransposeecon]);

                                        [Drvtwotwoecon,nDrvtwotwoecon,sDrvtwotwoecon] = lmivar(3, -1*sRecon);
                                        %
                                        %likewise R has no imaginary component
                                        %
                                        PRAecon = newlmi;

                                        %Structured singular value issue starts here
                                        lmiterm([1 0 0 0],Routrfctor) %
                                        lmiterm([1 1 1 Drvoneoneecon],1,1) %
                                        lmiterm([-1 1 1 Zecon],Aehat.',1,'s')

                                        lmiterm([1 1 2 Drvonetwoecon],1,1) %
                                        lmiterm([-1 1 2 Zecon],1,Behat) %

                                        lmiterm([1 2 2 Drvtwotwoecon],1,1) % 

                                        %Positive definitedness of
                                        %parahermitian issue

                                        lmiterm([2 1 1 V],Ad.',1,'s') %
                                        lmiterm([-2 1 2 Eecon],1,1) %
                                        lmiterm([2 1 2 V],1,Bd) %
                                        lmiterm([-2 2 2 Recon],1,1) %

                                        %positive definiteness
                                        %condition on the matrix
                                        %Rsprecon
                                        %
                                        lmiterm([-3 1 1 Recon],1,1) %for full Rsprecon > 0
                                        


                                        %
                                        lmisysecon = getlmis;

                                        %Determining the LMI feasibility solutino and assigning it to a matrix
                                        %storage
                                        [tminecon,xfeasecon] = feasp(lmisysecon,[0,200,0,25,0],tminlimit);
                                        if ((isempty(xfeasecon) == 0) && (tminecon < 0))
                                        %
                                        %Piecewise frequency extraction of the decision variables
                                        %
                                        Sopteconrsltdummy = dec2mat(lmisysecon,xfeasecon,Eecon);
                                        Ropteconrsltdummy = dec2mat(lmisysecon,xfeasecon,Recon);
                                        numecondecvr = decnbr(lmisysecon);

                                        %
                                        LargeZecondummy = dec2mat(lmisysecon,xfeasecon,Zecon);
                                        Vdummy = dec2mat(lmisysecon,xfeasecon,V);
                                        %--------------------------------------------------------------------------
                                        %Defining riccati equation parameters 
                                        %(Ad-Bd*inv(R)*S')'*Y1 + Y1*(Ad-Bd*inv(R)*S') - Y1*Bd*inv(R)*Bd'*Y1 + P - S*inv(R)*S'
                                        % where P is given as a zero matrix in this case from the related
                                        % optimisation.
                                        %--------------------------------------------------------------------------
                                        %
                                                        Popteconrslt = zeros(Adrow); 
                                                        Sopteconrslt = Sopteconrsltdummy;
                                                        Ropteconrslt = Ropteconrsltdummy;
                                                        %
                                                        Ahopteconrslt = Ad - Bd*inv(Ropteconrslt)*Sopteconrslt';
                                                        Bhopteconrslt = Bd*inv(Ropteconrslt)*Bd';
                                                        Qhopteconrslt = Popteconrslt - Sopteconrslt*inv(Ropteconrslt)*Sopteconrslt';

                                                        %
                                                        %Obtaining the stabilising solution
                                                        %

                                                        %Form the Hamiltonian matrix
                                                        Hammatrixecon = [Ahopteconrslt -Bhopteconrslt;-Qhopteconrslt -Ahopteconrslt'];
                                                        %
                                                        %Obtaining the stabilising solution
                                                        [Y21,Y22,failecon,reig_minecon] = ric_schr(Hammatrixecon);
                                                        if failecon == 0
                                                        Y2 = Y22*inv(Y21);

                                                        %
                                                        %increment counter
                                                        spctralcntecon = 1;
                                                        %
                                                        %
                                                        else
                                                        disp('Problem with the determination of the riccati stabilising solution as per Economised parameterisation')
                                                        end
                                        else
                                        disp('Problem with non-existence of a feasible solution as per Economised parameterisation')
                                        end
                                        %
                                        %
                                            Gammausednonrecursiveecon = gammachoiceb;
                                            gammachoiceb = gammachoiceb + 0.1;
                                        %
                                        %
                                        end
                                        %
                                        %
                                        %
                                        if ((spctralcntecon == 1) && (isempty(xfeasecon) == 0))
                                        %
                                        %
                                        %Stop timer in the existence of a
                                        %solution 
                                        Econtimer = toc;
                                        %
                                        %
                                        Necon = N;
                                        tauecon = tau;
                                        %
                                        %
                                        %
                                        disp(' ')
                                        disp(' --------------------------------------------------------------------------- ')
                                        disp(' ')
                                        disp(['Non-recursive Econ opt Feasible with gammab = ',num2str(Gammausednonrecursiveecon), ' , \tau = ', num2str(tauecon), ' and N = ', int2str(Necon)])
                                        disp(' ')
                                        disp(' --------------------------------------------------------------------------- ')
                                        disp(' ')
                                        %
                                        %
                                        %
                                        else
                                        %
                                        %Stop timer in case no solution was
                                        %obtained
                                        %
                                        Econdummytimer= toc;
                                        %
                                        %
                                        %
                                        end

                                       
                %
                %
                %                
                
             
                %-----------------------------------------------------------------
                %
                % DISPLAY OF THE RESULTS IN A USEFUL FORM STARTS FROM HERE
                %
                %
                %-----------------------------------------------------------------
                            
                 if (nonrecurmchktauNbrk  == 0)
                             %
                            if ((spctralcntfull == 1) && (spctralcntecon == 1) && (isempty(xfeasfull) == 0) && (isempty(xfeasecon) == 0))
                                
                            %
                            %Displaying the number of decision variables
                            %involved in each optimisation
                            %
                            
                            disp('-------------------------------------------------------------------------------------- ')
                            disp(' ')
                            disp(['Full decision variables optimisation had ', int2str(numfulldecvr ),' decision variables.'])
                            disp(' ')
                                                        
                            disp(' ')
                            disp(['Economised decision variables optimisation had ', int2str(numecondecvr),' decision variables.'])
                            disp(' ')
                            disp('-------------------------------------------------------------------------------------- ')
                            
                            %
                            %Displaying the time taken to run each segment of the
                            %programss
                            %

                            disp(' ')
                            disp(' ')
                            disp(['Full decision variables optimisation took ', num2str(Fulltimer),' seconds to complete using N = ', int2str(Nfull),' and tau = ', num2str(taufull), ' .'])
                            disp(' ')

                            disp(' ')
                            disp(['Economised decision variables optimisation took ', num2str(Econtimer),' seconds to complete using N = ', int2str(Necon),' and tau = ', num2str(tauecon), ' .'])
                            disp(' ')
                            
                            disp('-------------------------------------------------------------------------------------- ')

                            %
                            %Displaying the values
                            %
                            disp(' ')
                            disp(' ')
                            disp(' The symmetrix matrix Pd for full parameterisation is ')
                            disp(' ')
                            seesys(Poptfullrslt,'%18.8f')
                            disp('-------------------------------------------------------------------------------------- ')
                                                      
                            %                           
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp(' The symmetrix matrix Z from full parameterisation is ')
                            disp(' ')
                            seesys(LargeZdummy,'%18.8f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' The symmetrix matrix Z from economised parameterization is ')
                            disp(' ')
                            seesys(LargeZecondummy,'%18.8f')
                            disp('-------------------------------------------------------------------------------------- ')

                            disp(' ')
                            disp(' ')
                            disp(' The symmetrix matrix Y from full parameterizatoin is ')
                            disp(' ')
                            seesys(Jdummy,'%18.8f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' The symmetrix matrix Y from economised parameterizatoin is ')
                            disp(' ')
                            seesys(Vdummy,'%18.8f')
                            disp('-------------------------------------------------------------------------------------- ')
                            %
                            %
                            %
                            %
                            %Compare the matrix R of the economised and the full version
                            %
                            disp(' ')
                            disp(' ')
                            disp('From full parameterization, R was given as ')
                            seesys(Roptfullrslt,'%18.8f')
                            disp(' ')

                            disp(' ')
                            disp('From economised parameterization, R was given as ')
                            seesys(Ropteconrslt,'%18.8f')
                            disp(' ')
                            disp('-------------------------------------------------------------------------------------- ')
                            %
                            %Compare the matrix S of the economised and the full version
                            %
                            disp(' ')
                            disp(' ')
                            disp('From full parameterization, S was given as ')
                            seesys(Soptfullrslt,'%18.8f')
                            disp(' ')

                            disp(' ')
                            disp('From economised parameterization, Sde was given as ')
                            seesys(Sopteconrslt,'%18.8f')
                            disp(' ')
                            disp('-------------------------------------------------------------------------------------- ')
                            %
                            %
                            %
                            %-------------------------------------------------------------------------%
                            %
                            % Showing the transfer functions and Plotting the SVD of the original parahermitian and the economised
                            % parahermitian to compare and contrast.
                            %
                            %-------------------------------------------------------------------------%
                            %
                            % The common eigenvalues of both the full D(s) and Economised D(s) is give
                            % below
                            %
                            %FOR THE FULL SPECTRAL FACTORISATION
                            %
                            %
                            %Displaying the stabilising solution
                            %

                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('--- Stabilising solution for the Riccati Equation for full parameterization ---')
                            disp(' ')
                            seesys(Y1,'%18.8f')
                            %
                            %

                            %
                            %Using Ad, Bd and the Phoptfullrslt and Shoptfullrslt and Rhoptfullrslt obtained from the optimisation problem
                            %construct then apply spectra factorisation to obtain the stable minimum
                            %phase transfer function D(s)_full
                            %
                            dummy1 = Soptfullrslt';
                            dummy2 = Bd'*Y1;
                            Fhat = -1*sqrtm(Roptfullrslt)*(dummy1 + dummy2);
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Ffullrslt is ')
                            seesys(Fhat,'%18.8f')
                            Cd = -1*Fhat;
                            %
                            Dd = sqrtm(Roptfullrslt);
                            %
                            %Statespace Representation of the Unit transfer function
                            %
                            Dfull = pck(Ad,Bd,Cd,Dd);
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('D(s) for full decision variable parahermitian is given as:')
                            disp(' ')
                            seesys(Dfull,'%18.8f') 
                            disp('-------------------------------------------------------------------------------------- ')

                            %-----------------------------------------

                            %FOR THE ECONOMISED SPECTRAL FACTORISATION
                            %
                            %
                            %
                            %Displaying the stabilising solution
                            %
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('--- Stabilising solution for the Riccati Equation for Economised Option ---')
                            disp(' ')
                            seesys(Y2,'%18.8f')
                            %
                            %
                            %
                            %Using Ad, Bd and the Popteconrslt and Sopteconrslt and Ropteconrslt obtained from the optimisation problem
                            %construct then apply spectra factorisation to obtain the stable minimum
                            %phase transfer function D(s)_econ.
                            %
                            dummy1 = Sopteconrslt';
                            dummy2 = Bd'*Y2;
                            Fhat = -1*sqrtm(Ropteconrslt)*(dummy1 + dummy2);
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Feconrslt is ')
                            seesys(Fhat,'%18.8f')
                            Cd = -1*Fhat;
                            %
                            Dd = sqrtm(Ropteconrslt);
                            %
                            %Statespace Representation of the Economised Unit transfer function
                            %
                            Decon = pck(Ad,Bd,Cd,Dd);
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('D(s) for Economised decision variable parahermitian is given as:')
                            disp(' ')
                            seesys(Decon,'%18.8f')
                            disp('-------------------------------------------------------------------------------------- ')
                            %
                            %
                            %
                            disp(' ')
                            disp('The Eigenvalues of both D(s)_econ and D(s)_full are:')
                            disp(' ')
                            Ad_eigenvalues = eig(Ad);
                            Ad_eigenvalues
                            %
                            %Generating the frequency response data for the different parahermitian
                            %functions (i.e. economised and original one)
                            %
                            %
                            Dfullfrqrsp = frsp(Dfull,omega);
                            Deconfrqrsp = frsp(Decon,omega);
                            disp(' ')
                            disp(' ')
                            disp('System information for full D(s) is ')
                            minfo(Dfull)
                            disp(' ')
                            disp('System information for economised D(s) is ')
                            minfo(Decon)
                            disp(' ')
                            disp('-------------------------------------------------------------------------------------- ')

                            %
                            %Plotting the different comparison graphs between economic and original
                            %parameterisation
                            %

                            [uDfull,sDfull,vDfull] = vsvd(Dfullfrqrsp);
                            [uDecon,sDecon,vDecon] = vsvd(Deconfrqrsp);
                             
                            %
                            %Plotting the different comparison graphs between economic and original
                            %parameterisation
                            %

                            figure
                            vplot('liv,lm',sDfull,sDecon, 'x:')
                            ylabel('Singular Values of  D(s)_{{}_{full}} :  -    &    D(s)_{{}_{econ}} :  x','FontSize',11);
                            xlabel('Frequency (rad/sec)','FontSize',11);
                            title(['For  {\itN} = ',num2str(N),'   and   \tau = ',num2str(tau)],'FontSize',11)
                            print -depsc compareds
                            print -djpeg compareds
                            
                            %
                            %Constructing the parahermitian matrix functions both the full
                            %parameterization and the economised parameterizations
                            %
                            Aparafull = [Ad zeros(Adrow);-Poptfullrslt -Ad'] ;
                            Bparafull = [Bd ; -Soptfullrslt];
                            Cparafull = [Soptfullrslt' Bd'];
                            Dparafull = Roptfullrslt;
                            %
                            Parafull = pck(Aparafull,Bparafull,Cparafull,Dparafull);


                            %For the economised
                            Aparaecon = [Ad zeros(Adrow);zeros(Adrow) -Ad'] ;
                            Bparaecon = [Bd ; -Sopteconrslt];
                            Cparaecon = [Sopteconrslt' Bd'];
                            Dparaecon = Ropteconrslt;
                            %
                            Paraecon = pck(Aparaecon,Bparaecon,Cparaecon,Dparaecon);

                            %
                            %
                            %
                            
                            disp('System information for Full Parahermitian is ')
                            minfo(Parafull)
                            disp(' ')
                            disp('System information for Economised Parahermitian is ')
                            minfo(Paraecon)
                            disp(' ')
                            disp('-------------------------------------------------------------------------------------- ')

                            %
                            %Comparing the pointwise H-infinity norm of the parahermitian matrix function
                            %obtained
                            %
                            
                            figure
                            Paraeconvarmat = frsp(Paraecon,omega);
                            Parafullvarmat = frsp(Parafull,omega);
                            [uParaecon,sParaecon,vParaecon] = vsvd(Paraeconvarmat);
                            [uParafull,sParafull,vParafull] = vsvd(Parafullvarmat);
                            vplot('liv,lm',sParafull,sParaecon,'x:')
                            ylabel('Singular Values of  \Gamma(s)_{{}_{full}} :  -    &    \Gamma(s)_{{}_{econ}} :  x','FontSize',11);
                            xlabel('Frequency (rad/sec)','FontSize',11);
                            title(['For  {\itN} = ',num2str(N),'   and   \tau = ',num2str(tau)],'FontSize',11)
                            print -depsc comparegammas
                            print -djpeg comparegammas

                            figure
                            nurmParaecon = vnorm(Paraeconvarmat);
                            nurmParafull = vnorm(Parafullvarmat);
                            vplot('liv,m',nurmParafull,'c-',nurmParaecon,'yx')
                            axis([10^-4 10^4 0 3000])
                            gammalegend = legend('$${}\|{}\Gamma(s)_{{}_{full}}{}\|_{_{\infty}}$$','$${}\|{}\Gamma(s)_{{}_{econ}}{}\|_{_{\infty}}$$',0);
                            set(gammalegend,'Interpreter','none')
                            set(gammalegend,'Interpreter','latex')
                            whitebg([1 1 1])
                            ylabel('Pointwise Frequency plot of  {}||{}\Gamma(s)_{{}_{full}}{}||{}_{\infty}    and    {}||{}\Gamma(s)_{{}_{econ}}{}||{}_{\infty}','FontSize',11);
                            xlabel('Frequency (rad/sec)','FontSize',11);
                            title(['For  {\itN} = ',num2str(N),'   and   \tau = ',num2str(tau)],'FontSize',11)
                            print -depsc comparegammanorms
                            print -djpeg comparegammanorms
                            
                            disp(' ')
                            disp('The Eigenvalues of Parafull')
                            disp(' ')
                            Aparafull_eigenvalues = eig(Aparafull);
                            Aparafull_eigenvalues
                            
                            disp(' ')
                            disp('The Eigenvalues of Paraecon')
                            disp(' ')
                            Aparaecon_eigenvalues = eig(Aparaecon);
                            Aparaecon_eigenvalues
                            
                            %
                            %Calculate and compare the norms of the
                            %parahermitian matrix functions
                            %
                            [dataparaecon,rowpointparaecon,indvparaecon,errparaecon] = vunpck(nurmParaecon);
                            [dataparafull,rowpointparafull,indvparafull,errparafull] = vunpck(nurmParafull);
                            hinfnormdiffpara = abs(dataparafull - dataparaecon);
                            hinfnormdiffparavarmat = vpck(hinfnormdiffpara,indvparaecon);
                            
                            figure
                            vplot('liv,m',hinfnormdiffparavarmat, 'g-')
                            axis([10^-4 10^4 0 10])
                            gammalegend = legend('$$ \; | \; {} \; {}||{}\Gamma(s)_{{}_{full}}{}||{}_{\infty}  -  {}||{}\Gamma(s)_{{}_{econ}}{}||{}_{\infty} \; {} \; {} | \; $$',1);
                            set(gammalegend,'Interpreter','none')
                            set(gammalegend,'Interpreter','latex')
                            whitebg([1 1 1])
                            ylabel('Pointwise Frequency plot of   |  {}||{}\Gamma(s)_{{}_{full}}{}||_{_{\infty}}    -    ||{}\Gamma(s)_{{}_{econ}}{}||_{_{\infty}}{}  |','FontSize',11);
                            xlabel('Frequency (rad/sec)','FontSize',11);
                            title(['For  {\itN} = ',num2str(N),'   and   \tau = ',num2str(tau)],'FontSize',11)
                            print -depsc gammanormdiff
                            print -djpeg gammanormdiff
                            
                            %
                            %Approximate the given system matrix or constant
                            %matrix and then output the approximated form
                            %
                            Paraeconapprox = roundupsys(Paraecon,3);
                            Parafullapprox = roundupsys(Parafull,3);
                            %
                            Dfullapprox = roundupsys(Dfull,3);
                            Deconapprox = roundupsys(Decon,3);
                            %
                            Poptfullrsltapprox = roundupconst(Poptfullrslt,2);
                            Soptfullrsltapprox = roundupconst(Soptfullrslt,2);
                            Roptfullrsltapprox = roundupconst(Roptfullrslt,2);
                            %
                            Sopteconrsltapprox = roundupconst(Sopteconrslt,2); 
                            Ropteconrsltapprox = roundupconst(Ropteconrslt,2);
                            
                            [testa,testb,testc,testd] = unpck(Dfull);
                             testaapprox = roundupconst(testa,3);
                             testbapprox = roundupconst(testb,3);
                             testcapprox = roundupconst(testc,3);
                             testdapprox = roundupconst(testd,3);
                             
                            [testecona,testeconb,testeconc,testecond] = unpck(Decon);
                             testeconaapprox = roundupconst(testecona,3);
                             testeconbapprox = roundupconst(testeconb,3);
                             testeconcapprox = roundupconst(testeconc,3);
                             testecondapprox = roundupconst(testecond,3);
                            %
                            %
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original A_full is given as:')
                            disp(' ')
                            seesys(testa,'%19.10f')
                            disp(' ')
                            disp('Approximated A_full is given as:')
                            disp(' ')
                            seesys(testaapprox,'%15.3f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original B_full is given as:')
                            disp(' ')
                            seesys(testb,'%19.10f')
                            disp(' ')
                            disp('Approximated B_full is given as:')
                            disp(' ')
                            seesys(testbapprox,'%15.3f')
                            
                             disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original C_full is given as:')
                            disp(' ')
                            seesys(testc,'%19.10f')
                            disp(' ')
                            disp('Approximated C_full is given as:')
                            disp(' ')
                            seesys(testcapprox,'%15.3f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original D_full is given as:')
                            disp(' ')
                            seesys(testd,'%19.10f')
                            disp(' ')
                            disp('Approximated D_full is given as:')
                            disp(' ')
                            seesys(testdapprox,'%15.3f')
                            %
                            %
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original D(s)full is given as:')
                            disp(' ')
                            seesys(Dfull,'%19.10f')
                            disp(' ')
                            disp('Approximated D(s)full is given as:')
                            disp(' ')
                            seesys(Dfullapprox,'%15.3f')
                            
                            %
                            %
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original A_econ is given as:')
                            disp(' ')
                            seesys(testecona,'%19.10f')
                            disp(' ')
                            disp('Approximated A_econ is given as:')
                            disp(' ')
                            seesys(testeconaapprox,'%15.3f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original B_econ is given as:')
                            disp(' ')
                            seesys(testeconb,'%19.10f')
                            disp(' ')
                            disp('Approximated B_econ is given as:')
                            disp(' ')
                            seesys(testeconbapprox,'%15.3f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original C_econ is given as:')
                            disp(' ')
                            seesys(testeconc,'%19.10f')
                            disp(' ')
                            disp('Approximated C_econ is given as:')
                            disp(' ')
                            seesys(testeconcapprox,'%15.3f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original D_econ is given as:')
                            disp(' ')
                            seesys(testecond,'%19.10f')
                            disp(' ')
                            disp('Approximated D_econ is given as:')
                            disp(' ')
                            seesys(testecondapprox,'%15.3f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original D(s)econ is given as:')
                            disp(' ')
                            seesys(Decon,'%18.8f')
                            disp(' ')
                            disp('Approximated D(s)econ is given as:')
                            disp(' ')
                            seesys(Deconapprox,'%15.3f')
                            
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Gammafull is given as:')
                            disp(' ')
                            seesys(Parafull,'%18.8f')
                            disp(' ')
                            disp('Approximated Gammafull is given as:')
                            disp(' ')
                            seesys(Parafullapprox,'%15.3f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Gammaecon is given as:')
                            disp(' ')
                            seesys(Paraecon,'%18.8f')
                            disp(' ')
                            disp('Approximated Gammaecon is given as:')
                            disp(' ')
                            seesys(Paraeconapprox,'%15.3f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Pdfull is given as:')
                            disp(' ')
                            seesys(Poptfullrslt,'%19.9f')
                            disp(' ')
                            disp('Approximated Pdfull is given as:')
                            disp(' ')
                            seesys(Poptfullrsltapprox,'%16.2f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Sdfull is given as:')
                            disp(' ')
                            seesys(Soptfullrslt,'%19.9f')
                            disp(' ')
                            disp('Approximated Sdfull is given as:')
                            disp(' ')
                            seesys(Soptfullrsltapprox,'%16.2f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Sdecon is given as:')
                            disp(' ')
                            seesys(Sopteconrslt,'%19.9f')
                            disp(' ')
                            disp('Approximated Sdecon is given as:')
                            disp(' ')
                            seesys(Sopteconrsltapprox,'%16.2f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Rdfull is given as:')
                            disp(' ')
                            seesys(Roptfullrslt,'%18.9f')
                            disp(' ')
                            disp('Approximated Rdfull is given as:')
                            disp(' ')
                            seesys(Roptfullrsltapprox,'%15.2f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Recon is given as:')
                            disp(' ')
                            seesys(Ropteconrslt,'%18.9f')
                            disp(' ')
                            disp('Approximated Recon is given as:')
                            disp(' ')
                            seesys(Ropteconrsltapprox,'%15.2f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Ad*1000 is given as:')
                            disp(' ')
                            Adscaled = Ad*(10^2);
                            seesys(Adscaled,'%18.10f')
                            disp(' ')
                            disp('Approximated Ad*1000 is given as:')
                            disp(' ')
                            Adapproxscaled = roundupconst(Adscaled,2);
                            seesys(Adapproxscaled,'%16.2f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Bd*1000 is given as:')
                            disp(' ')
                            Bdscaled = Bd*(10^3);
                            seesys(Bdscaled,'%18.10f')
                            disp(' ')
                            disp('Approximated Bd*1000 is given as:')
                            disp(' ')
                            Bdapproxscaled = roundupconst(Bdscaled,2);
                            seesys(Bdapproxscaled,'%16.2f')
                            disp('-------------------------------------------------------------------------------------- ')
                            
                            [Aprint,Bprint,Cprint,Dprint] = unpck(Msys);
                            
                            disp(' ')
                            disp('M(s) state matrix A_g is given as: ')
                            disp(' ')
                            seesys(Aprint,'%15.2f')
                            
                            disp(' ')
                            disp('M(s) input matrix B_g is given as: ')
                            disp(' ')
                            seesys(Bprint,'%15.2f')
                            
                            disp(' ')
                            disp('M(s) Output matrix C_g is given as: ')
                            disp(' ')
                            seesys(Cprint,'%15.2f')
                            
                            disp(' ')
                            disp('M(s) matrix D_g is given as: ')
                            disp(' ')
                            seesys(Dprint,'%15.2f')
                            disp('-------------------------------------------------------------------------------------- ')
                            disp(' ')
                            disp(['The row of A1 is : ' int2str(tvra1)])
                            disp(' ')

                            disp(' ')
                            disp(['The row of A2 is : ' int2str(tvra2)])
                            disp(' ')

                            [tvrad,tvcad] = size(Ad);
                            disp(' ')
                            disp(['The row of Ad (i.e. row of A1 and row of A2) is : ' int2str(tvrad)])
                            disp(' ')
                            disp('-------------------------------------------------------------------------------------- ')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Rdfull is given as:')
                            disp(' ')
                            seesys(Roptfullrslt,'%18.9f')
                            disp(' ')
                            disp('Original Recon is given as:')
                            disp(' ')
                            seesys(Ropteconrslt,'%18.9f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Approximated Rdfull is given as:')
                            disp(' ')
                            seesys(Roptfullrsltapprox,'%15.2f')
                            disp(' ')
                            disp('Approximated Recon is given as:')
                            disp(' ')
                            seesys(Ropteconrsltapprox,'%15.2f')
                            disp('-------------------------------------------------------------------------------------- ')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original Sdfull is given as:')
                            disp(' ')
                            seesys(Soptfullrslt,'%19.9f')
                            disp(' ')
                            disp('Original Sdecon is given as:')
                            disp(' ')
                            seesys(Sopteconrslt,'%19.9f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Approximated Sdfull is given as:')
                            disp(' ')
                            seesys(Soptfullrsltapprox,'%16.2f')
                            disp(' ')
                            disp('Approximated Sdecon is given as:')
                            disp(' ')
                            seesys(Sopteconrsltapprox,'%16.2f')
                            disp('-------------------------------------------------------------------------------------- ')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original C_full is given as:')
                            disp(' ')
                            seesys(testc,'%19.10f')
                            disp(' ')
                            disp('Original C_econ is given as:')
                            disp(' ')
                            seesys(testeconc,'%19.10f')
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Approximated C_full is given as:')
                            disp(' ')
                            seesys(testcapprox,'%15.3f')
                            disp(' ')
                            disp('Approximated C_econ is given as:')
                            disp(' ')
                            seesys(testeconcapprox,'%15.3f')
                            disp('-------------------------------------------------------------------------------------- ')
                            
                            
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Original D_full is given as:')
                            disp(' ')
                            seesys(testd,'%19.10f')
                            disp(' ')
                            disp('Original D_econ is given as:')
                            disp(' ')
                            seesys(testecond,'%19.10f')
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp('Approximated D_full is given as:')
                            disp(' ')
                            seesys(testdapprox,'%15.3f')
                            disp(' ')
                            disp('Approximated D_econ is given as:')
                            disp(' ')
                            seesys(testecondapprox,'%15.3f')
                           
                            %Determining and showing the structure of the
                            %scaling transfer function matrix obtained
                            
                            %Obtain LTI model of Full D(s)
                            [astruct,bstruct,cstruct,dstruct] = unpck(Dfull);
                            Dfullstructdummy = ss(astruct,bstruct,cstruct,dstruct);
                            
                            %Obtain frequency response at of Full D(s) specific frequency points
                            Dfullstructone = freqresp(Dfullstructdummy,omega(1));
                            Dfullstructtwo = freqresp(Dfullstructdummy,omega(2));
                            Dfullstructthree = freqresp(Dfullstructdummy,omega(3));
                            Dfullstructfour = freqresp(Dfullstructdummy,omega(4));
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_full at w = ', num2str(omega(1)),' is given as: '])
                            disp(' ')
                            seesys(Dfullstructone,'%15.3f')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_full at w = ', num2str(omega(2)),' is given as: '])
                            disp(' ')
                            seesys(Dfullstructtwo,'%15.3f')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_full at w = ', num2str(omega(3)),' is given as: '])
                            disp(' ')
                            seesys(Dfullstructthree,'%15.3f')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_full at w = ', num2str(omega(4)),' is given as: '])
                            disp(' ')
                            seesys(Dfullstructfour,'%15.3f')
                            disp(' ')
                            
                            %Obtain LTI model of Economised D(s)
                            [astruct,bstruct,cstruct,dstruct] = unpck(Decon);
                            Deconstructdummy = ss(astruct,bstruct,cstruct,dstruct);
                            
                            %Obtain frequency response at of Economised D(s) specific frequency points
                            Deconstructone = freqresp(Deconstructdummy,omega(1));
                            Deconstructtwo = freqresp(Deconstructdummy,omega(2));
                            Deconstructthree = freqresp(Deconstructdummy,omega(3));
                            Deconstructfour = freqresp(Deconstructdummy,omega(4));
                            
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_econ at w = ', num2str(omega(1)),' is given as: '])
                            disp(' ')
                            seesys(Deconstructone,'%15.3f')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_econ at w = ', num2str(omega(2)),' is given as: '])
                            disp(' ')
                            seesys(Deconstructtwo,'%15.3f')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_econ at w = ', num2str(omega(3)),' is given as: '])
                            disp(' ')
                            seesys(Deconstructthree,'%15.3f')
                            disp(' ')
                            disp(' ')
                            disp(['The Structure of D_econ at w = ', num2str(omega(4)),' is given as: '])
                            disp(' ')
                            seesys(Deconstructfour,'%15.3f')
                            disp(' ')
                            
                            
                            %Use the frequency response data of Dfull and
                            %Decon to extract them and see their frequency
                            %dependent structure
                            
                            [Dfulldata,Dfullrowpoint,Dfullindv,Dfullerr] = vunpck(Dfullfrqrsp);
                            [Decondata,Deconrowpoint,Deconindv,Deconerr] = vunpck(Deconfrqrsp);
                            disp(' ')
                            disp(' ')
                            disp(' ')
                            %
                            %
                            %see the structure of Dfull
                            disp(' Structure of Dfull over the whole range of frequency is')
                            seesys(Dfulldata,'%15.3f')
                            disp(' ')
                            disp(' ')
                            %
                            %
                            %
                            %see the structure of Decon
                            disp(' Structure of Decon over the whole range of frequency is')
                            seesys(Decondata,'%15.3f')
                            disp(' ')
                            disp(' ')
                            %
                            %
                            save myvariables
                            save Dfull Dfull
                            save Decon Decon
                            %
                            disp(' ------------------------------------------------------- ')
                            disp(' ------------        DIARY ENDS       ------------------ ')
                            disp(' ------------------------------------------------------- ')
                            %
                            %
                            %
                            %
                            else
                            disp(' ')
                            disp(' ')
                            disp([' The range of values of Gamma and tau = ',num2str(tau),' and N = ',int2str(N),' did not provide the required solution'])
                            disp(' ')
                            disp(' ')
                            disp(' ------------------------------------------------------- ')
                            disp(' ------------        DIARY ENDS       ------------------ ')
                            disp(' ------------------------------------------------------- ')
                            end
                end
