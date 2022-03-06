% program to obtain the resolvent weights

clear
clc

%% input parameters
K1 = [-1            ,2            , 1];
K2 = [1            ,2           , 1];
K3 = [K1(1) + K2(1),K1(2) + K2(2), K1(3)];
Re = 200;
N  = 151;
Nsvd = 10;

[y, ~, ~, ~, U0] = channelMeanVel(Re, N);

X1 = 0:0.02:4*pi;
Y1 = -1:0.01:1;
Z1 = 0:0.1:1;

%% Create the set of 12 Fourier modes
triad1 = [K1;K2;K3];
triad2 = [-1*triad1(:,1), triad1(:,2), triad1(:,3)];
triad3 = [-1*triad1(:,1), -1*triad1(:,2), triad1(:,3)];
triad4 = [triad1(:,1), -1*triad1(:,2), triad1(:,3)];

triad = [triad1; triad2; triad3; triad4];

%% Get all interaction coefficients
% Obtain interaction coeffs with corresponding singvals and response modes

IntCoeff = zeros(1,6);
s = zeros(1,6);
psi = zeros(3*N,6);

[IntCoeff(1),s(1),psi(:,1)] = getIntCoeff(triad(3,:),triad(8 ,:),Re,N,Nsvd,U0,y);
[IntCoeff(2),s(2),psi(:,2)] = getIntCoeff(triad(3,:),triad(7 ,:),Re,N,Nsvd,U0,y);
[IntCoeff(3),s(3),psi(:,3)] = getIntCoeff(triad(1,:),triad(2 ,:),Re,N,Nsvd,U0,y);
[IntCoeff(4),s(4),psi(:,4)] = getIntCoeff(triad(6,:),triad(11 ,:),Re,N,Nsvd,U0,y);
[IntCoeff(5),s(5),psi(:,5)] = getIntCoeff(triad(6,:),triad(10,:),Re,N,Nsvd,U0,y);
[IntCoeff(6),s(6),psi(:,6)] = getIntCoeff(triad(4,:),triad(5 ,:),Re,N,Nsvd,U0,y);
disp(IntCoeff)

%% Create the system of equations and solve it
chi = 1e-7*[1,1,1,1,1,1];
re = 0;
im = 3;
if IntCoeff(1)*IntCoeff(3)>0 && IntCoeff(2)*IntCoeff(3)>0 && IntCoeff(2)*IntCoeff(1)>0 
    while (norm(chi)<1e-3)
        chi0 = (-10+20*rand(6,1) + 1j*(-10+20*rand(6,1)))';
        func = @(chi)weightEquations(chi,IntCoeff);
        options = optimoptions('fsolve','Algorithm','trust-region','Display','iter-detailed');
        chi = fsolve(func,chi0,options);
        disp(chi)
    end
else
    disp('constraints not satisfied')
end
%% create the optimization function and minimize it
% Problematic because it only optimizes over real values

% [x, D] = chebdif(N, 2);
% % D1 = D(:,:,1);
% D2 = D(:,:,2);
% ReStressMean = (D2*U0)/Re + 5;
% 
% f01 = getMeanf0(K1,Re,N,Nsvd,U0,y);
% f02 = getMeanf0(K2,Re,N,Nsvd,U0,y);
% f03 = getMeanf0(K3,Re,N,Nsvd,U0,y);
% 
% a = 0.99;
% fun = @(chi)(a*(abs(chi(3)*chi(2)'*IntCoeff(1) - chi(1)) +...
%                 abs(chi(3)*chi(1)'*IntCoeff(2) - chi(2)) +...
%                 abs(chi(1)*chi(2)'*IntCoeff(3) - chi(3)) +...
%                 abs(chi(6)*chi(5)'*IntCoeff(4) - chi(4)) +...
%                 abs(chi(6)*chi(4)'*IntCoeff(5) - chi(5)) +...
%                 abs(chi(4)*chi(5)'*IntCoeff(6) - chi(6)))^2 +...
%             (1-a)*(norm(ReStressMean - chi(1)*chi(1)'*f01(N+1:2*N) -...
%             chi(2)*chi(2)'*f02(N+1:2*N) - chi(3)*chi(3)'*f03(N+1:2*N))/norm(ReStressMean))^2);
% 
% [chi,fval] = fminunc(fun,real(chi0))

%% Reconstruct velocity
% A Fourier mode, say u1hat can be represented using resolvent modes as
% u1hat = chi(1)*s(1)*psi(1:N,1)
% total velocity field would be sum over all Fourier modes
% u = chi(1)*s1*psi1(1:N) + chi(2)*s2*psi2(1:N) +... chi(6)*s6*psi6(1:N)

[X, Y, Z] = meshgrid(X1,Y1,Z1);
uFinal = zeros(size(X));
vFinal = zeros(size(X));
wFinal = zeros(size(X));

% fpsiu = griddedInterpolant(flip(y),flip(psi(1:N,i)));
% plot(Y1,fpsiu())
fU0 = griddedInterpolant(flip(y),flip(U0));
uMean = fU0(Y);

for i=1:6
   kx = triad(i,1);
   kz = triad(i,2);
   fpsiu = griddedInterpolant(flip(y),flip(psi(1:N,i)));
   fpsiv = griddedInterpolant(flip(y),flip(psi(N+1:2*N,i)));
   fpsiw = griddedInterpolant(flip(y),flip(psi(2*N+1:3*N,i)));

   uFinal = uFinal + real(chi(i)*s(i)*fpsiu(Y).*exp(1i*kx*X + 1i*kz*Z));
   vFinal = vFinal + real(chi(i)*s(i)*fpsiv(Y).*exp(1i*kx*X + 1i*kz*Z));
   wFinal = wFinal + real(chi(i)*s(i)*fpsiw(Y).*exp(1i*kx*X + 1i*kz*Z));
end

uFinal = uFinal + uMean;

%% Plot the velocity field (in y)
figure(1)
plot(Y1,uFinal(:,5,10).*uFinal(:,5,10))

%% Plot the velocity field (in x y z)
figure(2)
party = 1:floor(length(Y1)/8);
data = squeeze(uFinal(party,:,1));
contourf(X1,party,data)
colormap("parula")
colorbar("eastoutside")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Possible sensible solutions

% (1,2,1)
% (-1,2,1)
% Re = 200
% N = 151
% -7.7328 - 0.3768i   7.7394 + 0.2028i  -0.0863 - 0.0065i   2.7303 - 7.2446i   7.6009 + 1.4712i   0.0454 - 0.0737i







