% channel resolvent in velocity-vorticity formulation 

clear
clc


%% input parameters
kx   = 1;
kz   = 6;
c    = 10;
om   = c*kx;
Re   = 2003;
N    = 201;
Nsvd = 16;

%% get mean velocity profile
[y, ~, ~, ~, U0] = channelMeanVel(Re, N);
y = y(2:N-1); % for no slip BCs

%% create operators
[L, C] = createOperators(kx, kz, N, Re, U0);

%% weight matrix - code taken from Barthel
k2    = kx^2 + kz^2;
Q     = generateWeightMatrix(N, k2);
sqrtQ = sqrtm(Q);

%% create resolvent operator and take its SVD
H = sqrtQ/(-1i*om*eye(2*(N-2)) - L)/sqrtQ;

[u, s, v] = svds(H, 20);
s = diag(s);
u = sqrtQ\u;
v = sqrtQ\v;

%% map to [u; v; w]
psi = C*u;
phi = C*v;

%% plot the results

figure
semilogy(s,'o');
ylabel('sigma'); xlabel('index');
title('Singular values (v-eta)');
% fname = sprintf('%d-%d-%d-%d-singular_values.png',Re,kx,kz,om);
% saveas(gcf,fname);

%% plot several psi (response) modes to see the trend
figure
Ni = N-2;

for i = 1:Nsvd
    subplot(Nsvd/2,2,i)
    plot(y, real(psi(1:Ni,i)), 'LineWidth', 2)
    ylabel('v'); xlabel('y');
end
sgtitle('Response modes') 
fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,om);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,fname);
