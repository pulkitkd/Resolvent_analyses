R=25;
L = [[-1/R, 1];[0, -2/R]];
omega = 0.01;
H = inv(1j*omega*eye(2) - L);
[U, S, V] = svd(H);
SS = zeros(10,1);
eigvals = eig(H);

% for t=0:119
%     Lt = [[-t/R, t];[0, -2*t/R]];
%     SS(t+1) = norm(exp(-Lt));
% end

% scatter(real(eigvals),imag(eigvals))
% 
for i = 1:10
    H = inv(1j*omega*eye(2) - L);
    [U, S, V] = svd(H);
    SS(i,1) = S(1,1);
    omega = omega + 0.1;
end

plot(SS)