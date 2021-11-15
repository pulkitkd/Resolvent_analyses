function vel = velocityField(u,s,chi,kx,kz,N)
% function to reconstruct velocity field from resolvent modes
% currently only using 1 resolvent mode

vel = s(1)*chi1*u*exp(1i*(kx*x + kz*z - omega*t));
end