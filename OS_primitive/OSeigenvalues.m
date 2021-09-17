clc
clear

% Orr-Somerfeld equation in primitive variable formulation 

%% input parameters

kx   = 1;
kz   = 1;
Re   = 5000;
N    = 150;

%% build operators

[y,~] = chebdif(N,2); 
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

[V,ee] = eig(L,M);
ee = -1i*diag(ee);
ee(isinf(real(ee))|isinf(imag(ee)))=-100-100i;
[~, ind] = sort(imag(ee),'descend');
ees = ee(ind);
Vs = V(:,ind);

%% plot the eigenvectors
subplot(1,2,1)
%branch A : 2,6,7,18
%branch P : 4,5,6,9,10
%branch S : 60+
%Squire(?): - 1,3,8,12,13
c1 = 6;
% c = ind(c1);
u = Vs(1:N , :);
v = Vs(N+1:2*N , :);
w = Vs(2*N+1:3*N , :);
p = Vs(3*N+1:4*N , :);

v1 = v(:,c1);
plot(y,real(v1),y,imag(v1),y,abs(v1),'linewidth',1.5);

%% plot the eigenvalues
subplot(1,2,2)
plot(ees,'.','markersize',12)
hold on
plot(real(ees(c1)),imag(ees(c1)),'ro')
axis([0 1 -1 0]);

