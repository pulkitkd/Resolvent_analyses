function [L, C] = buildOperators(kx, kz, N, Re, U0)

% Constructs LNS operator in velocity-vorticity formulation, as well as 
% outpout operator that maps to u, v, w. See Moarref et al., JFM 734, 2013
% for more details.

%% differentiation matrices
[~,DM] = chebdif(N,2); 
D1 = DM(:,:,1); % first derivative
D2 = DM(:,:,2); % second derivative

[~,D4] = cheb4c(N); % fourth derivatve; incorporates clamped boundary conditions

%% some definitions to use later
k2   = kx^2 + kz^2;
k4   = k2^2;
I    = eye(N-2);
Z    = zeros(N-2);
U0y  = D1*U0; % mean shear
U0yy = D2*U0; % second derivative of mean profile

%% incorporate homogenous BCs (Dirichilet)
D1   = D1(2:N-1,2:N-1);
D2   = D2(2:N-1,2:N-1);
U0   = U0(2:N-1);
U0y  = U0y(2:N-1);
U0yy = U0yy(2:N-1);

%% build Orr-Sommerfeld and Squire operators
U0   = diag(U0);
U0y  = diag(U0y);
U0yy = diag(U0yy);

Delta  = (D2 - k2*I); % Laplacian
Delta2 = (D4 - 2*k2*D2 + k4*I); % bilaplacian

L11 = Delta\(1/Re*Delta2 + 1i*kx*(U0yy - U0*Delta)); % O-S operator
L12 = Z;
L21 = -1i*kz*U0y; % coupling operator
L22 = 1/Re*Delta - 1i*kx*U0; % Sq operator
L   = [L11, L12; L21, L22];

%% also compute output operator that maps [v; eta] to [u; v; w]
Cu = [1i*kx*D1, -1i*kz*I]/k2;
Cv = [I, Z];
Cw = [1i*kz*D1, 1i*kx*I]/k2;
C  = [Cu; Cv; Cw];