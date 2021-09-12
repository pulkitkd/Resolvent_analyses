% TransGrowth
%
% Program that computes transient growth of the OSS for 
% flow cases Poiselle and Couette.
%
%
%
% INPUT 
%
% nosmod    = number of Orr-Sommerfeld modes
% R         = Reynolds number
% alp       = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2)  
% nosmod    = total number of modes for normal velocity
% iflag     = flag
%             iflag = 1: compute the maximum growth and 
%                        initial condition in time 
%                        interval [0,T]
%             iflag = 2: compute the initial disturbance 
%                        yielding maximum growth at time T
%
% OUTPUT 
% d         =   3D Orr-Sommerfeld matrix  
% M         =   energy matrix
%
    zi=sqrt(-1);

% input data

    iflow=input('Poiseuille (1) or Couette flow (2) ');
    nosmod=input('Enter N the number of OS modes: ');
    R=input('Enter the Reynolds number: ');
    alp=input('Enter alpha: ');
    beta=input('Enter beta: ');
    iflag=input(...
    '(1) Max growth in [Tmin,Tmax] (2) Max growth at T ');

    if iflag==1,
      Tmin=input('Enter Tmin: ');
      Tmax=input('Enter Tmax: ');
      T=[Tmin Tmax];
    else
      T=input('Enter T: ');
    end;

% generate Chebyshev differentiation matrices

    [D0,D1,D2,D4]=Dmat(nosmod);

% set up Orr-Sommerfeld matrices A and B 


    if iflow==1,
      [A,B,u]=pois(nosmod,alp,beta,R,D0,D1,D2,D4);
    else
      [A,B,u]=couet(nosmod,alp,beta,R,D0,D1,D2,D4);
    end;

% generate energy weight matrix

    ak2=alp^2+beta^2;
    M=energy(nosmod+1,nosmod+1,ak2);

% compute the Orr-Sommerfeld matrix (by inverting B)

    d=inv(B)*A;
    
% compute the optimal
    [flowin,flowot,gg]=optimal(d,T,M,ak2,iflag);
    figure(1); plot(gg(:,1),gg(:,2));grid on




  



