clear
clc

% channel resolvent in primitive variable formulation 

%% input parameters

K1 = [6, 6, 0.66];

kx   = K1(1);
kz   = K1(2);
c    = K1(3);
omega= c*kx;

Re   = 2003;
N    = 201;
Nsvd = 16; %(even number)

% x = linspace(0,1,N);
% z = linspace(0,1,N);
% t = zeros(1,);

%% get mean velocity profile
[y, ~, ~, ~, U0] = channelMeanVel(Re, N); %length(U0) = length(y)

%% get the resolvent matrix and take the SVD
[u, s, v] = getResolventSVD(kx, kz, omega, Re, N, Nsvd, U0);
[u, v] = channelSym(u, v, s, Nsvd);


