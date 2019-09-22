%
%clear screen, close all
%previous windows
%
close all
clear all
clear classes
clear java
clc
%
close all
clear all
clear classes
clear java
clc
%
%                            %
%load only the variables that you need
%
close all
clear all
clear classes
clear java
clc
%
diary wrksmu
disp(' ------------        NEW DIARY BEGINS    ----------------- ')
disp(' --------------------------------------------------------- ')

%
pack
%

load Dfull
load Decon
%
%
%
pack
%
%
%
%
%
%Publishing the results of the code execution
%publish('Latest','html')
%open html/Latest.html
%Getting necessary inputs
%
disp(' ')
disp(' ')
disp(' ')
%
%
plntinchkcnt = 0;
plntstatechkcnt = 0;
while plntinchkcnt < 1
numplantinput = input('What is the number of plant input ? ');
numplantoutput = input('What is the number of plant output  ? ');
if numplantinput == numplantoutput
plntinchkcnt = 1;
else
disp(' ')
disp('Please note that the plant must be square i.e. have the same number of inputs and outputs')
disp(' ')
end
end
%
while plntstatechkcnt < 1
numplantstate = input('What is the number of plant states  ? ');
if numplantstate >= 1
plntstatechkcnt = 1;
else
disp(' ')
disp('Please note that the plant is not allowed to have zero states')
disp(' ')
end
end
%
%Specifying the number of inputs and outputs to the controller
%
disp(' ')
disp(' ')
disp(' ')
disp(' ')
kinchkcnt = 0;
while kinchkcnt < 1
nmeaskinput = input('How many plant input do you want to come from the controller ? ');
ncontkoutput = input('How many plant output do you want to go into the controller ? ');
if ((nmeaskinput >= numplantoutput) && (ncontkoutput >= numplantinput))
disp(' ')
disp('Error: Controller inputs and output are both equal to or more than plant output and input respectively')
disp(' ')

elseif ((nmeaskinput >= numplantoutput) && (ncontkoutput <= 0))
disp(' ')
disp('Error: Controller inputs is equal to or more than plant output and controller ouput is less than or equal to zero')
disp(' ')

elseif ((nmeaskinput >= numplantoutput) || (nmeaskinput <= 0))
disp(' ')
disp('Error: Controller inputs is equal to or more than plant output or it is less than or equal to zero')
disp(' ')

elseif ((nmeaskinput <= 0) && (ncontkoutput <= 0))
disp(' ')
disp('Error: Both Controller inputs and output are less than or equal to zero ')
disp(' ')

elseif ((nmeaskinput <= 0) && (ncontkoutput >= numplantinput))
disp(' ')
disp('Error: Controller inputs is less than or equal to zero and outputs is greater than or equal to plant output ')
disp(' ')

elseif ((ncontkoutput >= numplantinput) || (ncontkoutput <= 0))
disp(' ')
disp('Error: Controller output greater than or equal to plant input or it is less than or equal to zero ')
disp(' ')

else
kinchkcnt = 1;
end
end

%
%Specifying the frequency response for the inputs
%
omega = logspace(-4,4,100);
%----------------------------------------------------------------------
% Initialising the plant G to be used to create the LFT M(s) by using a
% random 5 x 5 stable plant G note that for this program, i am using 
% an 8 state, 5 input and 5 output plant
%----------------------------------------------------------------------
rand('seed',0);
randn('seed',0);
Gdummy = rss(numplantstate,numplantoutput,numplantinput);
disp(' ')
disp('About to obtain a minimal realization for the plant')
disp(' ')
%Obtain the minimal realisation
dummy = minreal(Gdummy);
Gdummy = dummy;

%
%Extract the system matrix
%
[A,B,C,D] = ssdata(Gdummy);
%
%Create system matrix for the original plant and plot the singular values
%
%
Gdummysys = pck(A,B,C,D);
disp('The Original Plant has the minimal realization with ')
minfo(Gdummysys)
disp('and its packed form is')
seesys(Gdummysys,'%11.5f') 
disp(' ')
disp(' ')

%
figure
Gdummysysfrq = frsp(Gdummysys ,omega);
[ug,sGdummysysfrq,vg] = vsvd(Gdummysysfrq);
vplot('liv,lm',sGdummysysfrq);
ylabel('Singular values of G(s)','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -djpeg -append Gs
print -depsc Gs
%
%Obtain the continuous time state space model of the plant
%
Gctssm = ss(A,B,C,D);
%
%Synthesising the H-Infinity Controller and the LFT
%
[K,Tzw,vopt] = hinfsyn(Gctssm,nmeaskinput,ncontkoutput,0,10,0.001,2);
%
%Determining the frequency response of the controllers
%
gminlimit = 1.05*vopt;
[Ksub,Tzwsub,vsubopt] = hinfsyn(Gctssm,nmeaskinput,ncontkoutput,gminlimit,10,0.001,2);
disp(' ')
disp(' ')
disp(' Optimal controller synthesis ends ')
%
%Showing and Plotting the Controller
%
disp(' ')
disp(' ')
disp(' Suboptimal controller synthesis begins')
[Ak,Bk,Ck,Dk] = ssdata(K);
[Aksub,Bksub,Cksub,Dksub] = ssdata(Ksub);
%
% Pack into a system matrix
%
Ksys = pck(Ak,Bk,Ck,Dk);
disp(['The Optimal Controller provided an optimal value of ',num2str(vopt)])
minfo(Ksys)
disp('and its packed form is')
seesys(Ksys,'%11.5f') 
disp(' ')
disp(' ')
%
%
%
Ksubsys = pck(Aksub,Bksub,Cksub,Dksub);
disp(['The Sub-Optimal Controller provided an sub-optimal value of ',num2str(vsubopt)])
minfo(Ksubsys)
disp('and its packed form is')
seesys(Ksubsys,'%11.5f') 
disp(' ')
disp(' ')

%
%Obtain the frequency response
%
figure
Kfrq = frsp(Ksys,omega);
Ksubfrq = frsp(Ksubsys,omega);
[uk,sk,vk] = vsvd(Kfrq);
[uksub,sksub,vksub] = vsvd(Ksubfrq);
vplot('liv,lm',sk,'b--',sksub,'r--');
ylabel('Singular Value of  K(s) - (blue)   &   {K_{subopt}}(s) - (red)','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -depsc compareKs
print -djpeg compareKs
%
%Hence the Matrix M to be used in the D,K iteration is given as
%
M = Tzw;
Msub = Tzwsub;
%
%Extract the state matrix of the LFT
%
[Am,Bm,Cm,Dm] = ssdata(M);
[Amsub,Bmsub,Cmsub,Dmsub] = ssdata(Msub);


%
%Pack into a system matrix and scale by a multitude of 10
%

Msys = pck(Am,Bm,Cm,Dm);
Msys_old = Msys;
Msys = syscalarmult(Msys,1);
%Approximated and used in approximated form
Msys = roundupsys(Msys,2);

%Obtain the latest state-space matrices after approximation
[Am,Bm,Cm,Dm] = unpck(Msys);

%
%
disp(' ')
disp(' ------------------------------------------------------------------ ')
disp(' ')
%
%for Msys_old
[flagstabm,flagdectm] = stabdecttest(Msys_old);
%check the detectability and stabilizability of Msys_old(s)
if ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msys_old(s) is detectable and stabilizable')
disp(' ')
disp(' ----- ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msys_old(s) is detectable and but it is not stabilizable')
disp(' ')
disp(' ----- ')
elseif ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msys_old(s) is not detectable and but it is stabilizable')
disp(' ')
disp(' ----- ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msys_old(s) is neither detectable nor stabilizable')
disp(' ')
disp(' ----- ')
elseif ((strcmp(flagstabm,'Blank')== 1) || (strcmp(flagdectm,'Blank')== 1))
disp(' ')
disp('The detectability or stabilizability of the Msys_old(s) cannot be determined since it has no unstable modes ')
disp(' ')
disp(' ------ ')
end
%
%Testing the controllability and observability of Msys_old(s)  
%
[flagctrm,flagobsrm] = ctrobsrtest(Msys_old);
if ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msys_old(s) is observable and controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
elseif ((strcmp(flagctrm,'No')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msys_old(s) is observable and but it is not controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
elseif ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'No')== 1))
disp(' ')
disp('Msys_old(s) is not observable and but it is controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
else
disp(' ')
disp('Msys_old(s) is neither observable nor controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
end

%
%
%
%for Msys
[flagstabm,flagdectm] = stabdecttest(Msys);
%check the detectability and stabilizability of Msys(s)
if ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msys(s) is detectable and stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msys(s) is detectable and but it is not stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msys(s) is not detectable and but it is stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msys(s) is neither detectable nor stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'Blank')== 1) || (strcmp(flagdectm,'Blank')== 1))
disp(' ')
disp('The detectability or stabilizability of the Msys(s) cannot be determined since it has no unstable modes ')
disp(' ')
disp(' ------ ')
end
%
%Testing the controllability and observability of Msys(s)  
%
[flagctrm,flagobsrm] = ctrobsrtest(Msys);
if ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msys(s) is observable and controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
elseif ((strcmp(flagctrm,'No')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msys(s) is observable and but it is not controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
elseif ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'No')== 1))
disp(' ')
disp('Msys(s) is not observable and but it is controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
else
disp(' ')
disp('Msys(s) is neither observable nor controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
end

%
%For Msys and Msys_old
%
figure
Msys_oldfrq = frsp(Msys_old,omega);
Mfrq = frsp(Msys,omega);
[uMsys,sMsys,vMsys] = vsvd(Mfrq);
[uMsys_old,sMsys_old,vMsys_old] = vsvd(Msys_oldfrq);
vplot('liv,lm',sMsys_old,'b--',sMsys,'r--');
ylabel('Singular Value of  Msys_{old}(s) - (blue)   &   Msys(s) - (red)','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -depsc compareMsys_old
print -djpeg compareMsys_old
%
%
%

disp('The M lower lft matrix has ')
minfo(Msys)
disp('and its packed form is')
seesys(Msys,'%13.7f') 
disp(' ')
disp('and its eigenvalues are ')
eig(Am)
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
%Pack into a system matrix and scale by a multitude of 10
%

Msubsys = pck(Amsub,Bmsub,Cmsub,Dmsub);
Msubsys_old = Msubsys;
Msubsys = syscalarmult(Msubsys,1);
%Approximated and used in approximated form
Msubsys = roundupsys(Msubsys,2);

%Obtain the latest state-space matrices after approximation
[Amsub,Bmsub,Cmsub,Dmsub] = unpck(Msubsys);

%
%
%
%
%
disp(' ')
disp(' ------------------------------------------------------------------ ')
disp(' ')
%
%for Msubsys_old
[flagstabm,flagdectm] = stabdecttest(Msubsys_old);
%check the detectability and stabilizability of Msubsys_old(s)
if ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msubsys_old(s) is detectable and stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msubsys_old(s) is detectable and but it is not stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msubsys_old(s) is not detectable and but it is stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msubsys_old(s) is neither detectable nor stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'Blank')== 1) || (strcmp(flagdectm,'Blank')== 1))
disp(' ')
disp('The detectability or stabilizability of the Msubsys_old(s) cannot be determined since it has no unstable modes ')
disp(' ')
disp(' ------ ')
end
%
%Testing the controllability and observability of Msubsys_old(s)  
%
[flagctrm,flagobsrm] = ctrobsrtest(Msubsys_old);
if ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msubsys_old(s) is observable and controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
elseif ((strcmp(flagctrm,'No')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msubsys_old(s) is observable and but it is not controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
elseif ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'No')== 1))
disp(' ')
disp('Msubsys_old(s) is not observable and but it is controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
else
disp(' ')
disp('Msubsys_old(s) is neither observable nor controllable')
disp(' ')
disp(' +++++++++++++++++++++++++++++ ')
end

%
%
%
%for Msubsys
[flagstabm,flagdectm] = stabdecttest(Msubsys);
%check the detectability and stabilizability of Msubsys(s)
if ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msubsys(s) is detectable and stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'Yes')== 1))
disp(' ')
disp('Msubsys(s) is detectable and but it is not stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'Yes')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msubsys(s) is not detectable and but it is stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'No')== 1) && (strcmp(flagdectm,'No')== 1))
disp(' ')
disp('Msubsys(s) is neither detectable nor stabilizable')
disp(' ')
disp(' ------ ')
elseif ((strcmp(flagstabm,'Blank')== 1) || (strcmp(flagdectm,'Blank')== 1))
disp(' ')
disp('The detectability or stabilizability of the Msubsys(s) cannot be determined since it has no unstable modes ')
disp(' ')
disp(' ------ ')
end
%
%Testing the controllability and observability of Msubsys(s)  
%
[flagctrm,flagobsrm] = ctrobsrtest(Msubsys);
if ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msubsys(s) is observable and controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
elseif ((strcmp(flagctrm,'No')== 1) && (strcmp(flagobsrm,'Yes')== 1))
disp(' ')
disp('Msubsys(s) is observable and but it is not controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
elseif ((strcmp(flagctrm,'Yes')== 1) && (strcmp(flagobsrm,'No')== 1))
disp(' ')
disp('Msubsys(s) is not observable and but it is controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
else
disp(' ')
disp('Msubsys(s) is neither observable nor controllable')
disp(' ')
disp(' ------------------------------------------------------------------ ')
end

%
%For Msubsys and Msubsys_old
%
figure
Msubsys_oldfrq = frsp(Msubsys_old,omega);
Msubfrq = frsp(Msubsys,omega);
[uMsubsys,sMsubsys,vMsubsys] = vsvd(Msubfrq);
[uMsubsys_old,sMsubsys_old,vMsubsys_old] = vsvd(Msubsys_oldfrq);
vplot('liv,lm',sMsubsys_old,'b--',sMsubsys,'r--');
ylabel('Singular Value of  Msubsys_{old}(s) - (blue)   &   Msubsys(s) - (red)','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -depsc compareMsubsys_old
print -djpeg compareMsubsys_old
%
%
%
disp('The Msub lower lft matrix has ')
minfo(Msubsys)
disp('and its packed form is')
seesys(Msubsys,'%13.7f') 
disp(' ')
disp('and its eigenvalues are ')
eig(Amsub)
disp(' ')
%
%
%Approximate to a predetermined decimal places the state-space matrix of Msub
%
%
Msubdpsys = roundupsys(Msubsys,2);
%
%
%
%specifically done to crosscheck the approximation done using system matrix
disp('The Msubdp lower lft matrix is the approx of Msub and it has ')
minfo(Msubdpsys)
disp('and its packed form is')
seesys(Msubdpsys,'%13.2f') 
disp(' ')

[aMsubdpsys,bMsubdpsys,cMsubdpsys,dMsubdpsys] = unpck(Msubsys);
aMsubdpsys = roundupconst(aMsubdpsys,2);
disp('aMsubdpsys ')
seesys(aMsubdpsys,'%13.2f')
bMsubdpsys = roundupconst(bMsubdpsys,2);
disp('bMsubdpsys ')
seesys(bMsubdpsys,'%13.2f')
cMsubdpsys = roundupconst(cMsubdpsys,2);
disp('cMsubdpsys ')
seesys(cMsubdpsys,'%13.2f')
dMsubdpsys = roundupconst(dMsubdpsys,2);
disp('dMsubdpsys ')
seesys(dMsubdpsys,'%13.2f')
disp(' --------------------------------------------------------- ')
%
%Testing diary function
%

%
%Obtain the frequency response
%
figure
Mfrq = frsp(Msys,omega);
Msubfrq = frsp(Msubsys,omega);
[u,s,v] = vsvd(Mfrq);
[usub,ssub,vsub] = vsvd(Msubfrq);
vplot('liv,lm',s,'g--',ssub,'y--');
ylabel('Singular Value of  Final M(s) - (green)   &   Final {M_{subopt}}(s) - (yellow)','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -depsc compareMs
print -djpeg compareMs

%
%
%Choice on whether the user wants to use an optimal controller or a
%suboptimal controller
%
stringcnt = 0;
while stringcnt < 1
disp(' ')
disp(' ')
optstring1 = input('Do you want an Optimal controller to construct M(s)? Y/N: ','s');
if ((isempty(optstring1) == 0) && ((strncmpi('y', optstring1, 1) == 1) ||(strncmpi('n', optstring1, 1) == 1)))
    if (strncmpi('y', optstring1, 1) == 1)
    stringcnt = 1;
    disp(' ')
    disp(' You have choosen to use an optimal controller')
    disp(' ')
    else
    disp(' ')
    disp(' ')
    optstring2 = input('Do you instead want a Sub-optimal controller to construct M(s)? Y/N: ','s');
        if ((isempty(optstring1) == 0) && ((strncmpi('y', optstring2, 1) == 1) ||(strncmpi('n', optstring2, 1) == 1)))
            if (strncmpi('y', optstring2, 1) == 1)      
            disp(' ')
            disp(' You have choosen to use a  sub-optimal controller')
            disp(' ')
 %
 %Making the choice and switching from using an optimal controller to a
 %suboptimal controller
 %
            Mfrq = Msubfrq;
            Msys = Msubsys;
            stringcnt = 1;
            end
        else
            disp(' ')
            disp(' ')
            disp('Please answer the questions properly by either making a choice or using the correct alphabet')
            stringcnt = 0;
        end
    end
else
disp(' ')
disp(' ')
disp('Please answer the questions properly by either making a choice or using the correct alphabet')
stringcnt = 0;
disp(' ')
disp(' ')
end
end
%
%
pack
%
%
Mfrq = frsp(Msys,omega);

%
%
%%----------------------------------------------------------------------
%Obtaining the bounds for a default block structure for ease of
%illustration Mu synthesis command
%----------------------------------------------------------------------
% blk = [2 2;4 4];
 blk = [2 2;4 0]; %value of uncertainty assigned in the paper
% blk = [2 0;4 0]; %value of uncertainty assigned in the paper
[bounds,rowd] = mu(Mfrq,blk,'ltuC');

%
%Unpacking the bounds vector
%
[bnddat,rowpoint,indv,err] = vunpck(bounds);

%
% Extracting the upper and lower bound vectors for use later in the LMI
%
gammatestupper = bnddat(:,1);
gammatestuppermax = max(gammatestupper); 
disp(['The maximum upper bound is ', num2str(gammatestuppermax),' while the complete vectors of piecewise frequency upper bound is '])
seesys(gammatestupper','%10.5f')
disp(' ')
%
gammatestlower = bnddat(:,2);
disp('The complete vectors of piecewise frequency lower bound is ')
seesys(gammatestlower','%10.5f')
disp(' ')

%
%Determining and plotting the singular values
%using the obtained Dscale structures
%
dummydfull = minv(Dfull);
testdummy1 = mmult(Dfull,Msys,dummydfull);
testdummy1frq = frsp(testdummy1,omega);


%
dummydecon = minv(Decon);
testdummy2 = mmult(Decon,Msys,dummydecon);
testdummy2frq = frsp(testdummy2,omega);

%
%
%Plotting the bounds
%
figure
vplot('liv,lm', vnorm(testdummy1frq),'b--',bounds);
grid on
ylabel('\bar{\sigma}{}(D_{full}MD^{-1}{}_{full})','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -depsc DfullMDfullinv
print -djpeg DfullMDfullinv

%
figure
vplot('liv,lm',vnorm(testdummy2frq),'b--',bounds);
grid on
ylabel('\bar{\sigma}{}(D_{econ}MD^{-1}{}_{econ})','FontSize',11);
xlabel('Frequency (rad/sec)','FontSize',11);
print -depsc DeconMDeconinv
print -djpeg DeconMDeconinv

disp('-------------------------------------------------------------------------------------- ')
disp(' ------------                      NEW DIARY ENDS                   ------------------ ')