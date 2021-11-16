clear
clc

%% input parameters
K1 = [6, 6, 0.667];
K2 = [1, 6, 0.667];
K3 = [K1(1) + K2(1), K1(2) + K2(2), K1(3)];
Re   = 187;
N    = 201;

% Nsvd should be an even number
% Nsvd less than 10 can result in warnings stating inaccuracy in svds()
Nsvd = 10;

% The mean velocity remains the same for all interaction coeff.
[~, ~, ~, ~, U0] = channelMeanVel(Re, N);

triad1 = [K1;K2;K3];
triad2 = [-1*triad1(:,1), triad1(:,2), triad1(:,3)];
triad3 = [-1*triad1(:,1), -1*triad1(:,2), triad1(:,3)];
triad4 = [triad1(:,1), -1*triad1(:,2), triad1(:,3)];

triadFull = [triad1; triad2; triad3; triad4];

%% Interaction coefficient

IntCoeff(1) = getIntCoeff(triadFull(3,:),triadFull(8,:),Re,N,Nsvd,U0);
IntCoeff(2) = getIntCoeff(triadFull(3,:),triadFull(7,:),Re,N,Nsvd,U0);
IntCoeff(3) = getIntCoeff(triadFull(1,:),triadFull(2,:),Re,N,Nsvd,U0);
IntCoeff(4) = getIntCoeff(triadFull(6,:),triadFull(11,:),Re,N,Nsvd,U0);
IntCoeff(5) = getIntCoeff(triadFull(6,:),triadFull(10,:),Re,N,Nsvd,U0);
IntCoeff(6) = getIntCoeff(triadFull(4,:),triadFull(5,:),Re,N,Nsvd,U0);
IntCoeff(7) = getIntCoeff(triadFull(9,:),triadFull(2,:),Re,N,Nsvd,U0);
IntCoeff(8) = getIntCoeff(triadFull(9,:),triadFull(1,:),Re,N,Nsvd,U0);
IntCoeff(9) = getIntCoeff(triadFull(7,:),triadFull(8,:),Re,N,Nsvd,U0);
IntCoeff(10) = getIntCoeff(triadFull(12,:),triadFull(5,:),Re,N,Nsvd,U0);
IntCoeff(11) = getIntCoeff(triadFull(12,:),triadFull(4,:),Re,N,Nsvd,U0);
IntCoeff(12) = getIntCoeff(triadFull(10,:),triadFull(11,:),Re,N,Nsvd,U0);
