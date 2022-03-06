clear
clc
% channel resolvent in primitive variables
%% input parameters
kx   = 1;
kz   = 0;
c    = 1;
omega= c*kx;
Re   = 200;
N    = 101;
Nsvd = 8; %(even number)

% Reynolds-Tiederman mean velocity
[y, ~, ~, ~, U0] = channelMeanVel(Re, N); %length(U0) = length(y)
% U0 = 1-y.^2;
% plot(y,U0,'linewidth',2)
yc = getCritLayer(y,U0,c);

%% get the resolvent matrix and take the SVD
[u, s, v] = getResolventSVD(kx, kz, omega, Re, N, Nsvd, U0);

%% enforce channel symmetry
% [u, v] = channelSym(u, v, s, Nsvd);

%% plot the results
field = 3;
plotResolventModes(y, u, s, v, N, Nsvd, field, yc, Re, kx, kz, omega);

%% Validity of first order approximation
% [W, iW] = weight_matrix(N);
% H = u*diag(s)*v';
% H1 = u(:,1)*s(1)*v(:,1)';
% disp("Validity of first order approximation =")
% norm(H1-H)

%% mean reynolds stress from mean velocity
% mean(u'.grad(u')) = -u . grad u + grad p + laplacian u / Re
% figure
% [x, D] = chebdif(N, 2);
% D1 = D(:,:,1);
% D2 = D(:,:,2);
% ReStress = (D2*U0)/Re + 1;
% plot(y,ReStress)

%% Additional comments
% To verify orthogonality, take the inner product as follows
% W comes from the weight_matrix function as follows
% disp("Checking orthogonality")
% 
% W = W(1:3*N,1:3*N);
% testOrtho1 = u'*(W'*W)*u;
% norm(testOrtho1)
% 
% testOrtho2 = v'*(W'*W)*v;
% norm(testOrtho2)

% % Critical layer height
% yc = getCritLayer(y,U0,c);
% 
% field = 0;
% figure()
% plot(y, real(u(field*N+1:(field+1)*N,2)), 'LineWidth', 2)
% xticks([-yc,yc])
% yticks([])
% xlabel('y')
% ylabel('u')
% grid on
% ax.LineWidth = 1.5;


