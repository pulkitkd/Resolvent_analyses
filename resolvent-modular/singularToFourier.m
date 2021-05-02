function [fourierFields] = singularToFourier(sigmas, chis, psi, nsvd)
%singularToFourier - sum up singular modes to create the Fourier mode
%this function uses the singular response modes, singular values and mode
%weights to approximate the Fourier mode. 

%inputs: 
% sigmas- 1D array of n ordered singular values
% chi   - 1D array of n ordered modal weights
% psi   - 2D array with each column as an ordered response mode
% nsvd  - no of SVD modes to use in the approximation

N = length(psi)/4
fourierFields = zeros(4*N, 1);
u = zeros(N);
w = zeros(N);
v = zeros(N);
p = zeros(N);
if length(sigmas) == length(chis)
    disp('dimensions of sigma and chi match. Continue.')
else
    disp('error! sigma and chi dont have the same size. Taking default.')
    chis = ones(length(sigmas), 1);
end

for i = 1 : nsvd
    u = psi(1:N,i)'/max(psi(1:N,i));
    v = psi(N+1:2*N,i)'/max(psi(N+1:2*N,i));
    w = psi(2*N+1:3*N,i)'/max(psi(2*N+1:3*N,i));
    p = psi(3*N+1:4*N,i)'/max(psi(3*N+1:4*N,i));
    psi(:,i) = [u v w p]';
    fourierFields = fourierFields + sigmas(i)*chis(i)*psi(:,i);
end



end

