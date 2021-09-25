function [L, C] = createOperators(kx, kz, N, Re, U0)

%% differentiation matrices
% Note that for the second derivative, one should use D2 and not D1*D1 when
% cropped matrices are being used. 


[~,DM] = chebdif(N,2);
D1 = DM(:,:,1); %(N x N)
D2 = DM(:,:,2); %(N x N)

% fourth derivatve with clamped BCs
[~,D4] = cheb4c(N); %(N-2 x N-2)

Z = zeros(N-2);
I = eye(N-2);
ksq = kx^2 + kz^2;
U0y = D1*U0; 
U0yy = D2*U0;

D1 = D1(2:N-1, 2:N-1); % N-1 x N-1
D2 = D2(2:N-1, 2:N-1); % N-1 x N-1
U0 = U0(2:N-1); % N-1
U0y = U0y(2:N-1); % N-1
U0yy = U0yy(2:N-1); % N-1


U0 = diag(U0); % N-1 x N-1
U0y = diag(U0y); % N-1 x N-1
U0yy = diag(U0yy); % N-1 x N-1

% delta = (ksq*I - D2);

% Los = 1i*kx*U0*delta + 1i*kx*Uyy + (1/Re)*(D4 - 2*ksq*D2 + ksq^2*I);
% Lsq = 1i*kx*U0 + (1/Re)*delta;
% 
% L11 = Los;
% L12 = Z;
% L21 = 1i*kx*Uy;
% L22 = Lsq;
k2 = kx^2 + kz^2;
k4 = k2^2;
Delta  = (D2 - k2*I); % Laplacian
Delta2 = (D4 - 2*k2*D2 + k4*I); % bilaplacian

L11 = Delta\(1/Re*Delta2 + 1i*kx*(U0yy - U0*Delta)); % O-S operator
L12 = Z;
L21 = -1i*kz*U0y; % coupling operator
L22 = 1/Re*Delta - 1i*kx*U0; % Sq operator

L = [L11 L12; L21 L22];

cond(L)

%% also compute output operator that maps [v; eta] to [u; v; w]
Cu = [1i*kx*D1, -1i*kz*I]/k2;
Cv = [I, Z];
Cw = [1i*kz*D1, 1i*kx*I]/k2;
C  = [Cu; Cv; Cw];

% M = [delta Z; Z I];

end