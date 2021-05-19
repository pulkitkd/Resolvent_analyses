% Re = 7.5182E+05;
Re = 75000
k = 1;
n = 10;
omega = 0.8;
N = 128;
nsvd = 20;
mode = 1:3;
[r,su,ss,sv,U0,yP,UP,dU0,dr,H] = resolventSVD(Re,k,n,omega,N,nsvd);
[usu uss usv] = svds(H,nsvd);
uss = diag(uss);
uvmode1sq = usu(N+1:2*N,mode).*conj(usu(N+1:2*N,mode));
uvmode1sq = uvmode1sq./max(uvmode1sq);
svmode1sq = su(N+1:2*N,mode).*conj(su(N+1:2*N,mode));
svmode1sq = svmode1sq./max(svmode1sq);
vmode3sq = usu(N+1:2*N,3).*conj(usu(N+1:2*N,3));
vmode3sq = vmode3sq/max(vmode3sq);

clf

hold on
plot(r,svmode1sq(:,1),'-b')
plot(r,svmode1sq(:,2),'--b')
plot(r,svmode1sq(:,3),':b')
% plot(r,uvmode1sq(:,1),'-r')
% plot(r,uvmode1sq(:,2),'--r')
% plot(r,uvmode1sq(:,3),':r')
title('v velocity modes')
legend('scaled','','','unscaled')
image2 = gcf
filename = sprintf('%d-%d-%d-%d-%d-%d-v_velocity_modes.png',Re,k,n,omega,N,nsvd);
exportgraphics(image2,filename);
hold off

semilogy(ss(1:20), '*','linewidth',4)
title('Singular values of the resolvent')
xlabel('index')
ylabel('singular value')

% clf('reset')

% semilogy(ss, '-ob');
% %hold on
% %loglog(uss, '-or')
% title('Singular values')
% legend('scaled','unscaled')
% image1 = gcf;
% filename = sprintf('%d-%d-%d-%d-%d-%d-singular_values.png',Re,k,n,omega,N,nsvd);
% exportgraphics(image1,filename);
% hold off



