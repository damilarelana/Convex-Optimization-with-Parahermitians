                %
                %
                clear all
                close all
                clc
                %
                %
                %
                diary independentmu
                %
                %
                
                disp(' ------------------------------------------------- ')
                disp(' Diary of Independent Mu ')
                disp(' ------------------------------------------------- ')
                
                %
                %
                %
                
                Ag = [-19.48 -0.91 -3.15 7.94 -4.57 -4.64; -1.02 -10.59 -0.64 23.82 -13.71 -13.92; -3.24 -0.70 -30.99 15.42 -8.88 -9.01; 0 -0.15 -0.12 -5.29 -6.24 -7.00; 0 2.14 1.72 14.97 -24.21 -14.59; 0 2.70 2.17 10.21 -11.55 -39.76];
                Bg = [1.25 -0.38 -1.27 0 -0.02 0.76; -0.64 -0.19 0.98 0 0.62 0; 0.58 -1.38 -0.04 0.86 0.47 0; 0 0.03 0 0 -0.05 0; 0 -0.36 0 0 0.71 0; 0 -0.45 0 0 0.89 0];
                Cg = [-23.77 8.66 4.30 13.44 -7.73 -7.85; -2.74 0.24 0.19 17.37 -10.00 -10.15; 0 -9.23 -24.44 0 0 0; 0 -0.54 -10.85 11.91 -6.85 -6.96; -5.11 1.44 11.20 -2.70 1.55 1.58; -0.02 0.16 5.95 11.73 -6.75 -6.86];
                Dg = [4.14 8.55 0 9.24 0.61 -5.50; 0 -0.40 3.52 0 0.79 0.40; 0 -7.26 11.33 0 0 0; 3.18 -4.72 0 0 0.54 0; 0 -6.07 0 10.28 -0.12 0; 0 -0.27 -0.52 3.95 0.53 11.53];
                M = pck(Ag,Bg,Cg,Dg);
                
                %Obtain frequency response
                
                omega = logspace(-4,4,100);
                Mfrq = frsp(M,omega); 
                
                %----------------------------------------------------------------------
                %Obtaining the bounds for a default block structure for ease of
                %illustration Mu synthesis command
                %----------------------------------------------------------------------
                blk = [2 0;4 0]; %value of uncertainty assigned in the paper
                [bounds,rowd] = mu(Mfrq,blk,'ltuC');
   
                %
                %Unpacking the bounds vector
                %
                [bnddat,rowpoint,indv,err] = vunpck(bounds);

                %
                %Plotting the bounds
                %
                figure
                vplot('liv,lm',bounds);
                ylabel('Structured Singular Value of M(s)','FontSize',11);
                xlabel('Frequency (rad/sec)','FontSize',11);
                print -depsc Msmuforthesis
                print -djpeg Msmuforthesis

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