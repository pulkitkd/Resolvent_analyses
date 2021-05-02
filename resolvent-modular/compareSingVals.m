% get singular values for a set of modes
Re = 7.5182E+05;
k = 0.1;
n = 3;
omega = 0.1;
N = 96;
nsvd = 2;
itr = 100;
singvalsS = zeros(itr,1);
singvalsU = zeros(itr,1);

for i = 1:itr
    k = k + 0.01;
    [r,su,ss,sv,U0,yP,UP,dU0,dr,H] = resolventSVD(Re,k,n,omega,N,nsvd);
    singvalsS(i) = ss(1);
    [usu uss usv] = svds(H,nsvd);
    uss = diag(uss);
    singvalsU(i) = uss(1);
    disp(i);
end

loglog(singvalsS, '-ob')
hold on
loglog(singvalsU, '-or')
title('Singular values')
legend('scaled','unscaled')
xlabel('k')
ylabel('largest singular value')

image1 = gcf
filename = sprintf('%d-%d-%d-%d-%d-singular_values_comparison_with_k.png',Re,n,omega,N,nsvd);
exportgraphics(image1,filename);

hold off