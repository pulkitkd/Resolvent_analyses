% plot the eigenvalues of the Orr Sommerfeld equations

% inputs
Re = 10000;
Ny = 256;
kx = 1;
kz = 0;

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
U = 1. - y.*y;
Uy = D * U;
Uyy = D2 * U;

% create operators
M = [[ksq*Idn - D2, Z]; [Z, Idn]];
Los = 1j*kx*(diag(U))*(ksq*Idn - D2) + diag(1j*kx*Uyy) + Reinv*(ksq*ksq*Idn + D4 - 2*ksq*D2);
Lsq = 1j*kx*diag(U) + Reinv*(ksq*Idn - D2);
L = [[Los, Z];[(1j*kz)*diag(Uy), Lsq]];

% apply boundary conditions

% v = dv/dy = eta = 0 at walls

L(1, :) = 0.;
L(n, :) = 0.;
L(1, 1) = 1.;
L(n, n) = 1.;

L(n+1, :) = 0.;
L(2*n, :) = 0.;
L(n+1, n+1) = 1.;
L(2*n, 2*n) = 1.;

L(2, :) = 0.;
L(n-1, :) = 0.;
L(2, 1:n) = D(1, :);
L(n-1, 1:n) = D(n, :);

M(1, :) = 0.;
M(2, :) = 0.;

M(n, :) = 0.;
M(n-1, :) = 0.;

M(n+1, :) = 0.;
M(2*n, :) = 0.;

% Evaluate the eigenvalues and plot

[V, E] = eig(L, M);
E = -1j*diag(E);
E = sort(E);
E = E(1:100);
scatter(real(E),imag(E))