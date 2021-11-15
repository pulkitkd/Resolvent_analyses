clear
clc

%% input parameters
K1 = [6, 6, 0.667];
K2 = [1, 6, 0.667];

kx   = [K1(1), K2(1)];
kz   = [K1(2), K2(2)];
c    = [K1(3), K2(3)];
Re   = 187;
N    = 201;
Nsvd = 16; %(even number)

omega= c.*kx;
[~,DM] = chebdif(N,2); 
D1 = DM(:,:,1); % first derivative - N X N
Z = zeros(N);

[W, iW] = weight_matrix(N);
W = W(1:3*N,1:3*N);
% test=u'*(W'*W)*u;

%% get mean velocity profile
[y, ~, ~, ~, U0] = channelMeanVel(Re, N); %length(U0) = length(y)

%% get the resolvent modes
[psi1, s1, phi1] = getResolventSVD(kx(1), kz(1), omega(1), Re, N, Nsvd, U0);
[psi1, phi1] = channelSym(psi1, phi1, s1, Nsvd);

[psi2, s2, phi2] = getResolventSVD(kx(2), kz(2), omega(2), Re, N, Nsvd, U0);
[psi2, phi2] = channelSym(psi2, phi2, s2, Nsvd);

[psi3, s3, phi3] = getResolventSVD(kx(1)+kx(2), kz(2), omega(2), Re, N, Nsvd, U0);
[psi3, phi3] = channelSym(psi3, phi3, s3, Nsvd);
% plotResolventModes(y, psi2, s2, phi2, N, Nsvd, 0);

%% compute the interaction coeff
% *TODO* - Ensure that we are picking the right even and odd modes!!!
% v should be odd, u and w even
%%
u1hat = psi1(1:N,2);
v1hat = psi1(N+1:2*N,2);
w1hat = psi1(2*N+1:3*N,2);

f3hat = phi3(:,1);

u2hat = psi2(1:N,2);
v2hat = psi2(N+1:2*N,2);
w2hat = psi2(2*N+1:3*N,2);

kx0 = kx(1)+kx(2);
kz0 = kz(1)+kz(2);

divpsipsi = 1i*kx0*eye(3*N)*[u1hat.*u2hat; u1hat.*v1hat; u1hat.*w1hat]+...
     [D1,Z,Z;Z,D1,Z;Z,Z,D1]*[v1hat.*u2hat; v1hat.*v2hat; v1hat.*w2hat]+...
     1i*kz0*eye(3*N)*[w1hat.*u2hat; w1hat.*v1hat; w1hat.*w1hat];

innerprod = f3hat'*(W'*W)*divpsipsi;
intcoeff = s1(1)*s2(1)*innerprod;