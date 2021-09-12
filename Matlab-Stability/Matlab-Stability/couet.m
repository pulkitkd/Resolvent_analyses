
function [A,B,u]=couet(nosmod,alp,beta,R,D0,D1,D2,D4);
%
% Function to create Orr-Sommerfeld matrices using Chebyshev 
% pseudospectral discretization for plane Couette flow 
% profile
%
% nosmod  = number of even or odd modes
% alp     = alpha
% R       = Reynolds number
% D0      = zero'th derivative matrix
% D1      = first derivative matrix
% D2      = second derivative matrix
% D4      = fourth derivative matrix

    zi=sqrt(-1);

% mean velocity

    ak2=alp^2+beta^2;
    Nos=nosmod+1;
    Nsq=nosmod+1;
    vec=(0:1:nosmod)';
    y=cos(pi*vec/nosmod);
    u=y;;
    du=ones(length(u),1);

%%%%%  Comments by Espen %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The below system gives an eigenvalue problem
% on the form  -omega*B + A=0, where the scalar 
% omega of course are eigenvalues (with corresponding 
% eigenvectors) of  the matrix inv(B)*A.
%%%%% End comments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% set up Orr-Sommerfeld matrix

    B11=D2-ak2*D0;
    A11=-(D4-2*ak2*D2+(ak2^2)*D0)/(zi*R);    % fourth derivative term
    A11=A11+alp*(u*ones(1,length(u))).*B11;   % add up the 
    er=-200*zi;
    A11=[er*[D0(1,:); D1(1,:)]; A11(3:Nos-2,:); ... 
         er*[D1(Nos,:); D0(Nos,:)]];
    B11=[D0(1,:); D1(1,:); B11(3:Nos-2,:); ... 
         D1(Nos,:); D0(Nos,:)];

% set up Squire matrix and cross-term matrix

    A21=beta*(du*ones(1,length(u))).*D0;
    A22=alp*(u*ones(1,length(u))).*D0-(D2-ak2*D0)/(zi*R);
    B22=D0;
    A22=[er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];
    
% combine all the blocks 

% A is the matrix containing the LOS (A11) LSQ (A22) and the 
% coupling term (A21)
A=[A11 zeros(Nos,Nsq); A21 A22];
% B is the first block matrix, containing the second derivative
B=[B11 zeros(Nos,Nsq); zeros(Nsq,Nos) B22];




