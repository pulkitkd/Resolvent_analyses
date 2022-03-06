function [f0] = getMeanf0(K,Re,N,Nsvd,U0,y)
% A function to compute the interaction coefficient between K and -K

kx  = K(1);
kz  = K(2);
c   = K(3);

omega= c.*kx;

[~,DM] = chebdif(N,2); 
D1 = DM(:,:,1); % first derivative - N X N
Z = zeros(N);

%% get the resolvent modes
[psi1, s1, phi1] = getResolventSVD(kx, kz, omega, Re, N, Nsvd, U0);
[psi1, ~] = channelSym(psi1, phi1, s1, Nsvd);

[psi2, s2, phi2] = getResolventSVD(-kx, -kz, -omega, Re, N, Nsvd, U0);
[psi2, ~] = channelSym(psi2, phi2, s2, Nsvd);

%% compute the interaction coeff
% Ensure that we are picking the right even and odd modes
% v should be odd, u and w even
% when we plot all three fields, it is seen that 
% -first column of u and w is even
% -first column of v is odd
% ...conveniently ?

u1hat = psi1(1:N,1);
v1hat = psi1(N+1:2*N,1);
w1hat = psi1(2*N+1:3*N,1);

u2hat = psi2(1:N,1);
v2hat = psi2(N+1:2*N,1);
w2hat = psi2(2*N+1:3*N,1);

% Check that u1hat , w1hat are even and v1hat is odd
if (isevenodd(u1hat)~=0 || isevenodd(w1hat)~=0 || isevenodd(v1hat)~=1)
    disp('u v w are not appropriately even or odd')
    plot(y,real(u1hat))
    hold on
    plot(y,real(v1hat))
    plot(y,real(w1hat))
    hold off
end

% Check that u2hat , w2hat are even and v2hat is odd
if (isevenodd(u2hat)~=0 || isevenodd(w2hat)~=0 || isevenodd(v2hat)~=1)
    disp('u v w are not appropriately even or odd')
    plot(y,real(u2hat))
    hold on
    plot(y,real(v2hat))
    plot(y,real(w2hat))
    hold off
end

%% compute the interaction coefficient

% Fourier transformed divergence of (psi psi) - a vector quantity
divpsipsi = [D1,Z,Z; Z,D1,Z; Z,Z,D1]*[v1hat.*u2hat; v1hat.*v2hat; v1hat.*w2hat];

f0 = s1(1)*s2(1)*divpsipsi;

if max(abs(imag(f0))) > 1e-15
    disp("Error! Mean stress from fluctuations is not purely real.")
end

end