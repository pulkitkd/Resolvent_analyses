clear
clc

% channel resolvent in primitive variable formulation 

%% input parameters
kx   = 1;
kz   = 10;
c    = 0.5;
omega= c*kx;
Re   = 180;
N    = 201;
Nsvd = 6;

%% get mean velocity profile

[y, ~, ~, ~, U0] = channelMeanVel(Re, N); %length(U0) = length(y) = N

%% create the LNS operators

[L, M] = LNS_operators(kx, kz, N, Re, U0);

LHS = L - 1i*omega*M;
RHS = M;

%% apply BCs
% no slip boundary conditions
% no condition for pressure
LHS(1,:)     = 0;
LHS(N,:)     = 0;
LHS(N+1, :)  = 0;
LHS(2*N,:)   = 0;
LHS(2*N+1,:) = 0;
LHS(3*N,:)   = 0;


LHS(1,1)         = 1;
LHS(N,N)         = 1;
LHS(N+1,N+1)     = 1;
LHS(2*N,2*N)     = 1;
LHS(2*N+1,2*N+1) = 1;
LHS(3*N,3*N)     = 1;


RHS(1,:)     = 0;
RHS(N,:)     = 0;
RHS(N+1,:)   = 0;
RHS(2*N,:)   = 0;
RHS(2*N+1,:) = 0;
RHS(3*N,:)   = 0;

%% create the resolvent operator

H = LHS\RHS;

%% scale it

[W, iW] = weight_matrix(N);
H = W*H*iW;

%% perform the SVD of the resolvent

[u, s, v] = svds(H, 20);
s = diag(s);
u = W*u;

%% plot the results
% singular values
figure
semilogy(s,'o');
ylabel('sigma'); xlabel('index');
title('Singular values (primitive)');
fname = sprintf('%d-%d-%d-%d-singular_values.png',Re,kx,kz,omega);
saveas(gcf,fname);

% singular vectors
field = 1; % u = 0, v = 1, w = 2, p = 3
figure
for i = 1:Nsvd
    subplot(Nsvd,1,i)
    plot(y, real(u(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
    ylabel('v'); xlabel('y');
end

sgtitle('Response modes') 
% fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% saveas(gcf,fname);

