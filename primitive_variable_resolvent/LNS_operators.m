function [L, M] = LNS_operators(kx, kz, N, Re)

% Constructs LNS operators

%% differentiation matrices
[y,DM] = chebdif(N,2); 
D1 = DM(:,:,1); % first derivative - N X N
D2 = DM(:,:,2); % second derivative - N X N

%% mean velocity
U0 = 1 - y.^2;

%% some definitions to use later
I    = eye(N);
Z    = zeros(N);
U0y  = D1*U0; % mean shear (N)

%% create the M operator

M = 1i*[I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z]; % 4N X 4N

%% create the L operator

Ld = 1i*kx*diag(U0) + (kx^2*I - D2 + kz^2*I)/Re; %(N X N)

L1 = [Ld      , diag(U0y) , Z       , 1i*kx*I]; %u 
L2 = [Z       , Ld        , Z       , D1     ]; %v 
L3 = [Z       , Z         , Ld      , 1i*kz*I]; %w
L4 = [1i*kx*I , D1        , 1i*kz*I , Z      ]; %contin.
L  = [L1; L2; L3; L4];

end