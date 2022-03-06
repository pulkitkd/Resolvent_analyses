function vel = velocityField(u,s,chi,kx,kz,N)
% function to reconstruct velocity field from resolvent modes
% currently using 3 resolvent modes for each Fourier mode

vel = s(1)*chi1*u*exp(1i*(kx*x + kz*z - omega*t));
end