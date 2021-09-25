clear
clc

% channel resolvent in primitive variable formulation 

%% input parameters
kx   = 1;
kz   = 10;
c    = 0.5;
omega= c*kx;
Re   = 2003;
N    = 201;
Nsvd = 16; %(even number)

%% get mean velocity profile

[y, ~, ~, ~, U0] = channelMeanVel(Re, N); %length(U0) = length(y)
% U0 = 1-y.^2;

%% create the LNS operators

[L, M] = LNS_operators(kx, kz, N, Re, U0);

LHS = (L - 1i*omega*M);
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

% Check the quadrature integral by hand calculation

[W, iW] = weight_matrix(N);
Hs = W*H*iW;
Hs = Hs(1:3*N,1:3*N);

%% take the SVD of the scaled resolvent

[su, s, sv] = svds(Hs, 20);
s = diag(s);
u = iW(1:3*N,1:3*N)*su;
v = iW(1:3*N,1:3*N)*sv;

%% enforce channel symmetry

[u, v] = channelSym(u, v, s, Nsvd);

%% plot the results
% singular values
figure
semilogy(s,'o');
ylabel('sigma'); 
xlabel('index');
title('Singular values (primitive)');
grid on
grid minor
% fname = sprintf('%d-%d-%d-%d-singular_values.png',Re,kx,kz,omega);
% saveas(gcf,fname);

% singular response modes
field = 2; % u = 0, v = 1, w = 2, p = 3

figure
for i = 1:Nsvd
    subplot(Nsvd/2,2,i)
    plot(y, real(u(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
    ylabel('v');
    xlabel('y');
end

sgtitle('Response modes (primitive)')
fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,fname);


% % singular forcing modes
% field = 1; % u = 0, v = 1, w = 2, p = 3
% 
% figure
% for i = 1:Nsvd
%     subplot(Nsvd/2,2,i)
%     plot(y, real(v(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
%     ylabel('v');
%     xlabel('y');
% end
% 
% sgtitle('Forcing modes')
% % fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% % saveas(gcf,fname);

