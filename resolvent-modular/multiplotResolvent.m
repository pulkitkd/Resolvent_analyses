% get all the fourier fields approximations from a sum of singular modes
nt = 4;
fourierFields = zeros(4*N, nt);
[response,singVal] = getSingularModes();
for i = 1:nt
    fourierFields(:,i) = singularToFourier(singVal(:,i), ones(2,1), response(:,(i-1)*nsvd+1:i*nsvd), nsvd);
end

U(:,it,ix)=u*exp(1i*k*x(ix)+1i*n*theta(it))*exp_omegat;