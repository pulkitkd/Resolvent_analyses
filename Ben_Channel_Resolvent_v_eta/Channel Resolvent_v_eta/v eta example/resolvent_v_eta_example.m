clear
clc

% channel resolvent in velocity-vorticity formulation 

%% input parameters
kx   = 1;
kz   = 5;
c    = 1.0;
om   = c*kx;
Re   = 180;
N    = 201;
Nsvd = 6;

%% get mean velocity profile
[y, ~, ~, ~, U0] = channelMeanVel(Re, N);
y = y(2:N-1); % for no slip BCs

%% build LNS and output operators
[L, C] = buildOperators(kx, kz, N, Re, U0);

%% build weight matrix to enforce energy norm
k2    = kx^2 + kz^2;
Q     = generateWeightMatrix(N, k2);
sqrtQ = sqrtm(Q);
% invsqrtQ = inv(sqrtQ);
%% compute resolvent operator
I = eye(2*(N-2));
% H = (-1i*om*I - L)\I;
H = sqrtQ/(-1i*om*I - L)/sqrtQ;
% H = sqrtQ*(-1i*om*I - L)*invsqrtQ;

%% SVD
[u, s, v] = svds(H, 20);
s = diag(s);
u = sqrtQ\u;
v = sqrtQ\v;

%% map to [u; v; w]
psi = C*u;
phi = C*v;

%% enforce symmetry in across channel centerline
[psi, phi] = channelSym(psi, phi, s, Nsvd);

%% plot psi_1
% Ni = N-2;
% figure
% subplot(3,1,1)
% plot(y,real(psi(Ni+1:2*Ni,1)),'LineWidth',2)
% ylabel('u');xlabel('y');
% subplot(3,1,2)
% plot(y,real(psi(1:Ni,1)),'LineWidth',2)
% ylabel('v');xlabel('y');
% subplot(3,1,3)
% plot(y,real(psi(2*Ni+1:3*Ni,1)),'LineWidth',2)
% ylabel('w');xlabel('y');

%% plot several psi (response) modes to see the trend
figure
Ni = N-2;

for i = 1:Nsvd
    subplot(Nsvd,1,i)
    plot(y, real(psi(1:Ni,i)), 'LineWidth', 2)
    ylabel('v'); xlabel('y');
end
sgtitle('Response modes') 
fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,om);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,fname);

%% plot several phi (forcing) modes to see the trend

figure;
Ni = N-2;

for i = 1:Nsvd
    subplot(Nsvd,1,i);
    plot(y, real(phi(1:Ni,i)), 'LineWidth', 2);
    ylabel('fv'); xlabel('y');
end

sgtitle('Forcing modes');
fname = sprintf('%d-%d-%d-%d-forcing_modes.png',Re,kx,kz,om);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,fname);

%% plot singular values
figure
semilogy(s,'o');
ylabel('sigma'); xlabel('index');
title('Singular values');
fname = sprintf('%d-%d-%d-%d-singular_values.png',Re,kx,kz,om);
saveas(gcf,fname);