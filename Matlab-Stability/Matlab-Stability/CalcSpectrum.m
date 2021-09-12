%
% a simple routine for calculation of eigenvalues and eigenvectors
% for the Orr-Sommerfeld-Squire equations (OSS). There are two baseflows
% implemented: 1) Plane Poiselle 2) Couette
%
% Functions (from the book Henningson&Schmidt (HS)):
% Dmat: 
%                    generate coefficients of differentiation matrices
% pois, couet:
%                    sets up OSS matrices for the two flowtypes
% iord2:
%                    computes eigenvalues and eigenvector  sorted by increasing
%                    imaginary part.
%
% Variables:
% thiscase = a string. either 'poiselle' or 'couette'
% nosmod  = number of modes (collocation points)
% alp     = alpha
% beta    = beta 
% R       = Reynolds number
% D0      = zero'th derivative matrix  [ T_0(y) T_1(y) ...   T_N(y)]  pp 484 
% D1      = first derivative matrix 
% D2      = second derivative matrix
% D3      = third derivative matrix 
% D4      = fourth derivative matrix

%%%%%%%%%%% Main variables
thiscase='poiselle';
R=10000;
alp=1;
beta=1;
nosmod=200;                       % number of collocation points

%%%% the complex variable
zi=sqrt(-1);

%%%% compute chebyshev points
vec=(0:1:nosmod)';
y=cos(pi*vec/nosmod);   

%%%% generate Chebyshev differentiation matrices (coefficients)
[D0,D1,D2,D4]=Dmat(nosmod);

%%%% set up Orr-Sommerfeld-Squire matrices A and B
if strcmp(thiscase,'poiselle')
  ax_imagmin=-1;
  ax_realmin=0;
  ax_realmax=1;
  [A,B,u]=pois(nosmod,alp,beta,R,D0,D1,D2,D4);
elseif  strcmp(thiscase,'couette')
  ax_imagmin=-1;
  ax_realmin=-1;
  ax_realmax=1;
  [A,B,u]=couet(nosmod,alp,beta,R,D0,D1,D2,D4);
end
disp(['Calculating eigenvalues for ',thiscase,' flow']) 
%%%% compute the Orr-Sommerfeld matrix (by inverting B)
d=inv(B)*A;

%%%% computes eigenvalues and sort by decreasing imag
[a_vecs,c_vals]=iord2(d);       % a_vecs and c_vals defined pp 489 and 

%%%%% transform eigenvectors to physical ([v eta]'=sum_n=0^N a_n T_n(y))
eig_vecs=[D0 zeros(size(D0)) ; zeros(size(D0)) D0]*a_vecs;

v=eig_vecs(1:length(y),:);
eta=eig_vecs(length(y)+1:end,:);


%%%%% plotting
figure(1)
plot(real(c_vals),imag(c_vals),'*')
xlabel('\omega_r'); ylabel('\omega_i')
axis([ax_realmin ax_realmax ax_imagmin 0.1])
grid on
title(['eigenmodes of the OSS-operator for alpha= ',num2str(alp),', beta= '...
       ,num2str(beta),', Re = ',num2str(R)])

while 1
inp=input('Give approx value of eigenvalue a+ib to plot (0 or 999 to exit) ')
if inp==0 || inp==999 || real(inp)>ax_realmax || real(inp)<ax_realmin ...
       || imag(inp)<ax_imagmin
  break
else
  
  [trash,ind]=min(c_vals-inp);
  inp=c_vals(ind);
end
figure(1)
hold on
circle=plot(real(c_vals(ind)),imag(c_vals(ind)),'ro');
hold off
figure(2)  
subplot(2,1,1)
plot(y,abs(v(:,ind)),y,real(v(:,ind)),y,imag(v(:,ind)))
title(['v, eigenvalue ',num2str(c_vals(ind))])
subplot(2,1,2)
plot(y,abs(eta(:,ind)),y,real(eta(:,ind)),y,imag(eta(:,ind)))
title(['\eta, eigenvalue ',num2str(c_vals(ind))])
end