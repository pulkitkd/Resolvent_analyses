clear
clc

% channel resolvent in primitive variables
%% input parameters
kx   = -7;
kz   = -12;
c    = 0.667;
omega= c*kx;
Re   = 187;
N    = 201;
Nsvd = 16; %(even number)

%% get mean velocity profile
[y, ~, ~, ~, U0] = channelMeanVel(Re, N); %length(U0) = length(y)

%% get the resolvent matrix and take the SVD
[u, s, v] = getResolventSVD(kx, kz, omega, Re, N, Nsvd, U0);
% [W, iW] = weight_matrix(N);
% W = W(1:3*N,1:3*N);

%% enforce channel symmetry
[u, v] = channelSym(u, v, s, Nsvd);
% test=u'*(W'*W)*u;

%% plot the results
plotResolventModes(y, u, s, v, N, Nsvd, 0);

