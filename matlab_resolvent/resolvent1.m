% plot the eigenvalues of the Orr Sommerfeld equations

% inputs
Re = 1000;
Ny = 128;
kx = 1;
kz = 6;
omega = 1;

Reinv = 1./Re;
n = Ny+1;
ksq = kx*kx + kz*kz;
n = Ny+1;

% Chebyshev matrices
[D, y] = cheb(Ny);
D2 = D*D;
D4 = D2*D2;
Z = zeros(n, n);
Idn = eye(n);
Id2n = eye(2*n);
U = 1. - y.*y;
Uy = D * U;
Uyy = D2 * U;

% create operators
M = [[ksq*Idn - D2, Z]; [Z, Idn]];
Los = 1j*kx*(diag(U))*(ksq*Idn - D2) + (1j*kx).*diag(Uyy) + Reinv*(ksq*ksq*Idn + D4 - 2*ksq*D2);
Lsq = (1j*kx).*diag(U) + Reinv*(ksq*Idn - D2);
L = [[Los, Z];[(1j*kz)*diag(Uy), Lsq]];

LHS = -1j*omega*M + L;
% apply boundary conditions

% v = dv/dy = eta = 0 at walls

LHS(1, :) = 0.;
LHS(n, :) = 0.;
LHS(1, 1) = 1.;
LHS(n, n) = 1.;

LHS(n+1, :) = 0.;
LHS(2*n, :) = 0.;
LHS(n+1, n+1) = 1.;
LHS(2*n, 2*n) = 1.;

LHS(2, :) = 0.;
LHS(n-1, :) = 0.;
LHS(2, 1:n) = D(1, :);
LHS(n-1, 1:n) = D(n, :);

% M(1, :) = 0.;
% M(2, :) = 0.;
% 
% M(n, :) = 0.;
% M(n-1, :) = 0.;
% 
% M(n+1, :) = 0.;
% M(2*n, :) = 0.;

H = inv(LHS);

[U, S, V] = svd(H);

s = diag(S);

semilogy(s(1:20), '*','linewidth',4)
title('Singular values of the resolvent')
xlabel('index')
ylabel('singular value')

v1 = abs(U(1:n,1));
v2 = abs(U(1:n,2));
v3 = abs(U(1:n,3));
v4 = abs(U(1:n,4));
figure()
subplot(4,1,1);
plot(y, v1)
subplot(4,1,2); 
plot(y, v2)
subplot(4,1,3); 
plot(y, v3)
subplot(4,1,4); 
plot(y, v4)