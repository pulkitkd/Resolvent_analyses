% get singular values for a set of modes
%Re = 7.5182E+05;
Re = 75000
k = 1;
n = 10;
omega = 0.1;
N = 128;
nsvd = 5;
itr = 50;
singvalsS = zeros(itr,nsvd);
singvalsU = zeros(itr,nsvd);

for i = 1:itr
    omega = omega + 0.02;
    [r,su,ss,sv,U0,yP,UP,dU0,dr,H] = resolventSVD(Re,k,n,omega,N,nsvd);
    singvalsS(i,:) = ss;
    [usu, uss, usv] = svds(H,nsvd);
    uss = diag(uss);
    singvalsU(i,:) = uss;
    disp(i);
end

loglog(singvalsS,'LineWidth',2.0)
% hold on
% loglog(singvalsU, '-or')
title('Singular values')
%legend('scaled','unscaled')
xlabel('omega')
ylabel('largest singular value')

image1 = gcf
filename = sprintf('%d-%d-%d-%d-%d-singular_values_comparison_with_omega.png',Re,k,n,N,nsvd);
exportgraphics(image1,filename);

hold off