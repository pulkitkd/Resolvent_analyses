% p40.m - eigenvalues of Orr-Sommerfeld operator (compare p38.m)
R = 5772; clf, [ay,ax] = meshgrid([.56 .04],[.1 .5]);
N = 100;
%2nd- and 4th-order differentiation matrices:
% [D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
% S = diag([0; 1 ./(1-x(2:N).^2); 0]);
% D4 = (diag(1-x.^2)*D^4 - 8*diag(x)*D^3 - 12*D^2)*S;
% D4 = D4(2:N,2:N);

% % differentiation matrices
[x,DM] = chebdif(N,2); 
D1 = DM(:,:,1); % first derivative
D2 = DM(:,:,2); % second derivative

[~,D4] = cheb4c(N); % fourth derivative; incorporates clamped boundary conditions

D1   = D1(2:N-1,2:N-1);
D2   = D2(2:N-1,2:N-1);

% Orr-Sommerfeld operators A,B and generalized eigenvalues:
I = eye(N-2);
A = (D4-2*D2+I)/R - 2i*I - 1i*diag(1-x(2:N-1).^2)*(D2-I);
B = D2-I;
[V,ee] = eig(A,B);
ee = diag(ee);
i = N/20-1; subplot('position',[ax(i) ay(i) .38 .38])
plot(ee,'.','markersize',12)
% axis([-.8 .2 -1 0])
% title(['N = ' int2str(N) ' \lambda_{max} = ' ...
% num2str(max(real(ee)),'%16.12f')]), drawnow
