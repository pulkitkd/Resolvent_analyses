clear
clc

% channel resolvent in primitive variable formulation 

%% input parameters

kx   = 1;
kz   = 0;
Re   = 2000;
N    = 100;

%% build operators
[y,~] = chebdif(N,1); 
[L, M] = LNS_operators(kx, kz, N, Re);

%% no slip boundary conditions
% no condition for pressure
L(1,:) = 0;
L(N,:) = 0;
L(N+1, :) = 0;
L(2*N,:) = 0;
L(2*N+1,:) = 0;
L(3*N,:) = 0;


L(1,1) = 1;
L(N,N) = 1;
L(N+1,N+1) = 1;
L(2*N,2*N) = 1;
L(2*N+1,2*N+1) = 1;
L(3*N,3*N) = 1;


M(1,:) = 0;
M(N,:) = 0;
M(N+1,:) = 0;
M(2*N,:) = 0;
M(2*N+1,:) = 0;
M(3*N,:) = 0;

[V,ee] = eig(L, M);
ee = diag(ee);
ee = sort(ee);

%% plot the eigenvalues
subplot(1,2,1)
plot(ee,'.','markersize',12)
axis([0 1 -1 0]);

subplot(1,2,2)
%% plot the eigenvectors
c = 90;
u1 = imag(V(1:N,c));
u1 = u1/max(abs(u1));
v1 = imag(V(N+1:2*N,c));
v1 = v1/max(abs(v1));
% p1 = imag(V(3*N+1:4*N,c));
% p1 = p1/max(abs(p1));
 
plot(y,u1)
hold on
plot(y,v1)
hold on
% plot(y,p1)
hold off

