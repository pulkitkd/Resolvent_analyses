% plot the energy in first singular mode relative to the total energy
% over a range of streamwise and spanwise wavenumbers
Re = 7.5182E+05;
omega = 0.5;
n1 = 1;
n2 = 100;
k1 = 0.01;
k2 = 10;

ksteps = 100;
nsteps = 10;

dk = k2 / ksteps;
dn = n2 / nsteps;
n = n1;
k = k1;

for i = 1:nsteps
   for j = 1:ksteps
      [r,su,ss,sv,U0,yP,UP,dU0,dr,H] = resolventSVD(Re,k,n,omega,N,nsvd);
      sigma1 = ss(1)
      totalE = ss(1)/dot(ss,ss)
      E(i,j) = sigma1/totalE
      lambdaxp = lambdaxp + dlx
   end
   lambdaxp = lambdax1p
   lambdazp = lambdazp + dlz
end