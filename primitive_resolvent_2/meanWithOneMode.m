clear
clc

% Check that we would get real velocity field if we kept all modes
% rather than discarding half of them

% Check with Ben's algorithm - fmincon()

% Min the cost functional (J(x)) defined as the weighted sum of the f(x) (the
% residual) + g(x) - > a norm of the difference between the two Re stress div
% calculations - i) the fluctuations ii) the meam mom eqs due to
% Re-Tiederman velocity profile. Then we have
% J(x) = a f(x) + (1-a) g(x)

% Current calculation is essentially a = 1. If a ~= 1, then we give some
% weight to g(x)

% input parameters
K1 = [1 ,6 ,0.667];

kx = K1(1);
kz = K1(2);
omega = K1(3)*kx;

Re = 1000;
N  = 251;
Nsvd = 10; %(even)
[y, ~, ~, ~, U0] = channelMeanVel(Re, N);

chi = 6; %(weight)

%grid to plot velocity field
X1 = 0:0.01:6*pi;
Y1 = -1:0.01:1;
Z1 = 0:0.1:1;
[X, Y, Z] = meshgrid(X1,Y1,Z1);

uFinal = zeros(size(X));
vFinal = zeros(size(X));
wFinal = zeros(size(X));

fU0 = griddedInterpolant(flip(y),flip(U0));
umean = fU0(Y);

[psi, s, phi] = getResolventSVD(kx, kz, omega, Re, N, Nsvd, U0);
[psi, phi] = channelSym(psi, phi, s, Nsvd);

fpsiu = griddedInterpolant(flip(y), flip(psi(1:N,1)));
fpsiv = griddedInterpolant(flip(y), flip(psi(N+1:2*N,1)));
fpsiw = griddedInterpolant(flip(y), flip(psi(2*N+1:3*N,1)));

uFinal = uFinal + real(chi*s(1)*fpsiu(Y).*exp(1i*kx*X + 1i*kz*Z));
vFinal = vFinal + real(chi*s(1)*fpsiv(Y).*exp(1i*kx*X + 1i*kz*Z));
wFinal = wFinal + real(chi*s(1)*fpsiw(Y).*exp(1i*kx*X + 1i*kz*Z));

uFinal = uFinal + umean;


% Plot the velocity field (in x y z)
figure
halfy = 1:floor(length(Y1)/8);
data = squeeze(uFinal(halfy,:,1));
contourf(X1,halfy,data)

figure
data = squeeze(uFinal(:,:,1));
contourf(X1,Y1,data)