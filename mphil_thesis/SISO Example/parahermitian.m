close all
clear all
clc
%
tic
% Given the following real matrices

A = [-6 -11 -6;1 0 0;0 1 0];
B = [1;0;0];
format short g
disp(' ')
disp(' ')
disp(' ')
disp('---Original given matrices are---')
disp(' ')
Ahat = A;
A
disp(' ')
%
Bhat = B;
B
disp(' ')
%---------------------------------------------------------------------%
%                   PROVE OF FIRST PART OF LEMMA 2                   %
%---------------------------------------------------------------------%
%
%
%
disp('---------------------------------------------------------------------')
disp(' ')
disp(' ')
disp('   COMMENCEMENT OF THE EXAMPLE FOR THE FIRST PART OF LEMMA 2     ')
disp(' ')
disp(' ')
disp('---------------------------------------------------------------------')
%
%
%
%To show that A is hurwitz
%
A_eigenvalues = eig(A);
disp('A is hurwitz because eigenvalues of A are given as')
A_eigenvalues
disp(' ')
disp('and also to ensure that A has not eigenvalues on the jw axis')
%
%Assuming values for matrix C and D to ensure a real rational stable and
%proper unit transfer function
%
C = [9 63 114];
D = 1;
P = C'*C;
R_t = D'*D;
S_t = C'*D;
%
% Obtaining the transfer function
%
[numt,dent]= ss2tf(A,B,C,D);

%Display the positive frequency response function in transfer function form
T = tf(numt,dent);
disp(' ')
disp(' ')
disp(' ')
disp('The unit transfer functn is given as:')
disp(' ')

T
T_zpk = zpk(T);
T_zpk
%
%Obtaining the inverse of the unit transfer function 
%
Atinv = A - B*inv(D)*C;
Btinv = B*inv(D);
Ctinv = -1*inv(D)*C;
Dtinv = inv(D);
%
% Obtaining the transfer function
%
[numtinv,dentinv]= ss2tf(Atinv,Btinv,Ctinv,Dtinv);

%Display the positive frequency response function in transfer function form
Tinv = tf(numtinv,dentinv);
disp(' ')
disp(' ')
disp(' ')
disp('The inverse of the unit transfer functn is given as:')
disp(' ')

Tinv
Tinv_zpk = zpk(Tinv);
Tinv_zpk

%
%
%
Atinv_eigenvalues = eig(Atinv);
disp(' ')
disp(' ')
disp(' ')
disp('Atinv is hurwitz because eigenvalues of Atinv are given as')
Atinv_eigenvalues
disp(' ')
disp('Therefore T is a unit transfer function since its inverse is also stable in RH-infinity')
%
%constructing the original parahermitian function
%

disp(' ')
disp(' ')
disp('------- constructing original parahermitian matrix function --------')
disp(' ')
disp(' ')
disp('With the following ')
disp(' ')
P
disp(' ')
R_t
disp(' ')
S_t
disp(' ')
disp(' ')
disp(' ')
disp('then the original positive definite Parahermitian transfer functn given below is shown to have no zeros on the jw axis:')
disp(' ')
%
%Setting up the control version for the old gamma to ensure that one
%obtained directly is correct
%
[arowsold,acolsold] = size(A);
Pdummy = -1*P;
Adummy = -1*A';
Sdummy1 = -1*S_t;
Sdummy2 = S_t';
Bdummy = B';
A_gammactrl = [A zeros(arowsold); Pdummy Adummy];
B_gammactrl = [B;Sdummy1];
C_gammactrl = [Sdummy2 Bdummy];
D_gammactrl = R_t;
%
%Showing the transfer function representation of the new parahermitian
%matrix function
%
Gamma = T'*T;
Gamma
Gamma_zpk = zpk(Gamma);
Gamma_zpk
%
disp(' ')
disp('- Control Version of the Old Gamma - ')
disp(' ')
[numgammactrl,dengammactrl] = ss2tf(A_gammactrl,B_gammactrl,C_gammactrl,D_gammactrl);
Gammactrl = tf(numgammactrl,dengammactrl);
Gammactrl
Gammactrl_zpk = zpk(Gammactrl);
Gammactrl_zpk
disp(' ')
disp(' ')
%
%
%Constructing the system matrices of the original parahermitian matrix
%function
%
%
[arows,acols] = size(A);
Pdummy = -1*P;
Adummy = -1*A';
Sdummy1 = -1*S_t;
Sdummy2 = S_t';
Bdummy = B';
A_gamma = [A zeros(arows); Pdummy Adummy];
B_gamma = [B;Sdummy1];
C_gamma = [Sdummy2 Bdummy];
D_gamma = R_t;
%
%Determining the Eigenvalues of the original parahermitian function
%
disp(' ')
disp(' ')
disp(' ')
disp('--- Displaying Eigenvalues of the Original Parahermitian Matrix Function ---')
disp(' ')

A_gamma_eigenvalues = eigs(A_gamma);
A_gamma_eigenvalues
%
% Constructing the Hamiltonian Matrix for the original parahermitian matrix
% function
%
H11 = A - B*inv(R_t)*S_t';
H12 = -1*B*inv(R_t)*B'; 
H21 = -1*(P - S_t*inv(R_t)*S_t'); 
H22 = -1*(A - B*inv(R_t)*S_t')'; 
%
disp(' ')
disp(' ')
disp(' ')
disp('--- Displaying the original Hamiltonian Matrix ---')
%
H_original = [H11 H12;H21 H22];
H_original
%
disp(' ')
disp(' ')
disp(' ')
disp('--- Displaying Eigenvalues of the Original Hamiltonian Matrix ---')
disp(' ')
H_oeigenvalues = eigs(H_original);
H_oeigenvalues
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('------- Implementing the result of Lemma 2.1 --------')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
%
% Now suppose Phat = 0 so as to test the result of the paper 
%
[rowp,colp] = size(P);
Phat = zeros(rowp);
Pdummy = P - Phat;
%
%Solving the Lyapunov equation
%
A_dummy = A';
X_t = lyap(A_dummy,Pdummy);
disp(' ')
disp(' ')
disp(' ')
disp('The unique real symmetric lyapunov solution is given as')
X_t
%
%Solving for Q12_t
%
Shat_t = S_t + X_t*B; 
%
Q12 = Shat_t;
Q22 = R_t;
%
%
%constructing the new parahermitian function by constructing it through its
%system matrices
%
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' constructing new parahermitian matrix function and testing its positive definiteness ')
disp(' ')
%
%
[arowsnew,acolsnew] = size(A);
Pdummy = -1*Phat;
Adummy = -1*A';
Sdummy1 = -1*Shat_t;
Sdummy2 = Shat_t';
Bdummy = B';
A_gammanew = [A zeros(arowsnew); Pdummy Adummy];
B_gammanew = [B;Sdummy1];
C_gammanew = [Sdummy2 Bdummy];
D_gammanew = R_t;
%
%Showing the transfer function representation of the new parahermitian
%matrix function
%
[numgammanew,dengammanew] = ss2tf(A_gammanew,B_gammanew,C_gammanew,D_gammanew);
Gammanew = tf(numgammanew,dengammanew);
disp(' ')
disp(' ')
disp(' ')
disp('The New parahermitian matrix function is given as:')
disp(' ')

Gammanew
Gammanew_zpk = zpk(Gammanew);
Gammanew_zpk

%
%Determining the Eigenvalues of the new parahermitian function
%
disp(' ')
disp(' ')
disp(' ')
disp('--- Displaying Eigenvalues of the New Parahermitian Matrix Function ---')
disp(' ')

A_gammanew_eigenvalues = eig(A_gammanew);
A_gammanew_eigenvalues
%
% Constructing the Hamiltonian Matrix for the original parahermitian matrix
% function
%
Hnew11 = A - B*inv(R_t)*Shat_t';
Hnew12 = -1*B*inv(R_t)*B'; 
Hnew21 = -1*(Phat - Shat_t*inv(R_t)*Shat_t'); 
Hnew22 = -1*(A - B*inv(R_t)*Shat_t')'; 
%
disp(' ')
disp(' ')
disp(' ')
disp('--- Displaying The new Hamiltonian Matrix ---')
disp(' ')
H_new = [Hnew11 Hnew12;Hnew21 Hnew22];
H_new
%
disp(' ')
disp(' ')
disp(' ')
disp('--- Displaying Eigenvalues of the New Hamiltonian Matrix ---')
disp(' ')
H_neweigenvalues = eig(H_new);
H_neweigenvalues
disp(' ')
disp(' ')
disp(' ')
disp('Now since the new HAMILTONIAN MATRIX has no eigenvalue on jw-axis,')
disp('it means that the new parahermitian matrix is positive definite (corrolary 13.16).')
disp('therefore it is seen that there exists')

Q12
disp(' ')
Q22
disp('and ?(jw)_hat > 0 such that result of Lemma 2.1 holds.')

%
%
%---------------------------------------------------------------------%
%                   PROVE OF SECOND PART OF LEMMA 2                   %
%---------------------------------------------------------------------%
%
disp(' ')
disp('---------------------------------------------------------------------')
disp(' ')
disp('    COMMENCEMENT OF THE EXAMPLE FOR THE SECOND PART OF LEMMA 2      ')
disp(' ')
disp('---------------------------------------------------------------------')
disp(' ')


Shat = Q12;
R = Q22;
%
%
disp(' ')
disp(' ')
disp(' ')
disp(' Q12 and Q22 obtaining via an LMI feasibility of the parahermitian ')
disp(' ')

Q12 = Shat;
Q12
%
disp(' ')
%
Q22 = R;
Q22
%--------------------------------------------------------------------------
%Defining the constant matrices to use for the determination fo the
%riccati stabilization solution to the riccati equation
%(A-B*inv(R)*Shat')'*Y2 + Y2*(A-B*inv(R)*Shat') - Y2*B*inv(R)*B'*Y2 - Shat*inv(R)*Shat'
%--------------------------------------------------------------------------

Atest = A - B*inv(R)*Shat';
Btest = B*inv(R)*B';
Qtest = -1*Shat*inv(R)*Shat';
%Form the Hamiltonian matrix
Ham = [Atest -Btest;-Qtest -Atest'];
%
%Obtaining the stabilising solution
%
%[Y2,L,G] = care(Atest,Btest,Qtest);
[Yham1,Yham2,fail,reig_min] = ric_schr(Ham);
if fail == 0
    Y2 = Yham2*inv(Yham1);
else
    disp('Problem with the determination of the riccati stabilising solution')
end

disp(' ')
disp(' ')
disp(' ')
disp('--- Stabilising solution for the Riccati Equation ---')
disp(' ')

Y2

%
%Applying spectra factorisation
%Start program for the spectra factorisation of the parahermitian matrix
%
dummy1 = Shat';
dummy2 = B'*Y2;
Fhat = -1*sqrtm(R)*(dummy1 + dummy2);
Fhat
Chat = -1*Fhat;
%
Dhat = sqrtm(R);

%
%Statespace Representation of the Unit transfer function
%
[nummhat,denmhat]= ss2tf(Ahat,Bhat,Chat,Dhat);
disp(' ')
disp(' ')
disp(' ')
disp('The State space representation of the Unit transfer functn is given as:')
disp(' ')
Ahat
%
Bhat
%
Chat
%
Dhat


%Display the positive frequency response function in transfer function form
Mhat = tf(nummhat,denmhat);

disp(' ')
disp('Transfer function representation of the Unit transfer functn is given as:')
disp(' ')

Mhat
Mhat_zpk = zpk(Mhat);
Mhat_zpk

disp(' ')
disp('Eigenvalue of the unit transfer functions are:')
disp(' ')

Ahat_eigenvalues = eig(Ahat);
Ahat_eigenvalues

%---------------------------------------------------------
%for the inverse version of the unit transfer function.
%---------------------------------------------------------
disp(' ')
disp('The State space representation of the Inverse of the Unit transfer functn is given as:')
disp(' ')

MhatinvA = Ahat + Bhat*inv(sqrtm(R))*Fhat;
MhatinvA
%
MhatinvB = Bhat*inv(sqrtm(R));
MhatinvB
%
MhatinvC = inv(sqrtm(R))*Fhat;
MhatinvC
%
MhatinvD = inv(sqrtm(R));
MhatinvD

%
%Showing the transfer function of the Mhatinverse
%
[nummhatinv,denmhatinv]= ss2tf(MhatinvA,MhatinvB,MhatinvC,MhatinvD);

%Display the positive frequency response function in transfer function form
Mhatinv = tf(nummhatinv,denmhatinv);

disp(' ')
disp('Transfer function representation of the Inverse of the Unit transfer functn is given as:')
disp(' ')
Mhatinv
Mhatinv_zpk = zpk(Mhatinv);
Mhatinv_zpk
%
disp(' ')
disp('Eigenvalue of the Inverse of the unit transfer functions are:')
disp(' ')

MhatinvA_eigenvalues = eig(MhatinvA);
MhatinvA_eigenvalues

%--------------------------------------------------------------------------
%
%Constructing the last parahermitian matrix function from this feasible values
%of Q12 and Q22 in order to test the result of Lemma 2.1
%
%--------------------------------------------------------------------------
disp(' ')
disp(' ')
disp(' ')
disp('Constructing The Last Parahermitian Matrix Function and Testing Its Positive Definiteness -')
disp(' ')
%
%
[arowslast,acolslast] = size(A);
Pdummy = -1*Phat;
Adummy = -1*A';
Sdummy1 = -1*Q12;
Sdummy2 = Q12';
Bdummy = B';
A_gammalast = [A zeros(arowsnew); Pdummy Adummy];
B_gammalast = [B;Sdummy1];
C_gammalast = [Sdummy2 Bdummy];
D_gammalast = Q22;
%
%Showing the transfer function representation of the new parahermitian
%matrix function
%
[numgammalast,dengammalast] = ss2tf(A_gammalast,B_gammalast,C_gammalast,D_gammalast);
Gammalast = tf(numgammalast,dengammalast);
disp(' ')
disp(' ')
disp(' ')
disp('The Last parahermitian matrix function is given as:')
disp(' ')

Gammalast
Gammalast_zpk = zpk(Gammalast);
Gammalast_zpk

%
%Determining the Eigenvalues of the last parahermitian function
%

disp(' ')
disp('--- Displaying Eigenvalues of the last Parahermitian Matrix Function ---')
disp(' ')

A_gammalasteigenvalues = eigs(A_gammalast);
A_gammalasteigenvalues
%
% Constructing the Hamiltonian Matrix for the last parahermitian matrix
% function
%
H11last = A - B*inv(Q22)*Q12';
H12last = -1*B*inv(Q22)*B'; 
H21last = -1*(Phat - Q12*inv(Q22)*Q12'); 
H22last = -1*(A - B*inv(Q22)*Q12')'; 
%
H_last = [H11last H12last;H21last H22last];
H_last
%
disp(' ')
disp(' ')
disp(' ')
disp('--- Displaying Eigenvalues of the New Hamiltonian Matrix ---')
disp(' ')
H_lasteigenvalues = eig(H_last);
H_lasteigenvalues
disp(' ')
disp(' ')
disp(' ')
disp('Now since this last HAMILTONIAN MATRIX has no eigenvalue on jw-axis,')
disp('it means that this last parahermitian matrix is positive definite (corrolary 13.16) as initially supposed.')
disp('therefore it is seen that there exists a unit transfer function T(jw) and ?(jw)_last > 0 such that result of Lemma 2.2 holds.')

%-------------------------------------------------------------------------%
%
% Plotting the SVD of the original parahermitian and the economised
% parahermitian to compare and contrast.
%
%-------------------------------------------------------------------------%
%
%Defining the frequency response ranges
%
omega = logspace(-4,4,100);
%
%Obtaining the System matrix Original Parahermitian (Gamma)
%
sysgamma = nd2sys(numgammactrl,dengammactrl);
%
%
%Obtaining the System matrix of Economised Parahermitian (Gammanew)
%
sysgammaeconmy = nd2sys(numgammanew,dengammanew);
%
%
%Obtaining the System matrix of the original unit transfer function
%
%
T_zpksys = nd2sys(numt,dent);
%
%Obtaining the System matrix of the Spectral Factorised unit transfer function
%
Mhat_zpksys = nd2sys(nummhat,denmhat);


%Generating the frequency response data for the different parahermitian
%functions (i.e. economised and original one)
%
%
sysgammafrqrsp = frsp(sysgamma,omega);
sysgammaeconmyfrqrsp = frsp(sysgammaeconmy,omega);
%
% Generating frequency response data for the original Unit transfer function and the Spectral Factorisation ones 
%
Mhat_zpkfrqrsp = frsp(Mhat_zpksys,omega);
T_zpkfrqrsp = frsp(T_zpksys,omega);
%
%Plotting the different comparison graphs between economic and original
%parameterisation
%

figure
vplot('liv,lm',sysgammafrqrsp,'bo-',sysgammaeconmyfrqrsp,'rx:')
ylabel('Magnitude')
xlabel('Frequency')
title('Plot of Original Parahermitian(Blue) & Economised Parahermitian(Red)')
%
%Plotting the different comparison graphs between original and spectral
%factorised unit transfer functions
%

figure 
vplot('liv,lm',T_zpkfrqrsp,'go-',Mhat_zpkfrqrsp,'yx:')
ylabel('Magnitude')
xlabel('Frequency')
title('Plot of Original Tf (Green) & Spectral Factor Tf(Yellow)')
toc