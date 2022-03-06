%% input parameters

clc
clear

kx = 1;
kz = 2;
Re = 4000;

Nsvd = 6; 

%% get mean velocity profile
N = 100;
[y, ~, ~, ~, U0] = channelMeanVel(Re, N); %length(U0) = length(y)
U0max = max(U0);
if kx == 0
    n = 30;
else
    n = ceil(kx)*ceil(U0max)+10*ceil(kx);
end
% figure
% plot(y,U0)

%% obtain the variation of singular values
singVals = zeros(n,4);
c = zeros(n,1);

for i = 1:n
    c(i) = i;
    [~, s, ~] = getResolventSVD(kx, kz, c(i), Re, N, Nsvd, U0);
    singVals(i,:) = s(1:4);
end

%% plot the results
figure
loglog(c/U0max, singVals(:,1),'-o')
hold on
loglog(c/U0max, singVals(:,2),'-o')
loglog(c/U0max, singVals(:,3),'-o')
loglog(c/U0max, singVals(:,4),'-o')
loglog(c/U0max, singVals(:,1)./singVals(:,3),'--')
grid on
hold off
xlabel('c/U_{0max}')
ylabel('largest singular value')
str = sprintf('Largest singVal with c for Re=%d kx=%.2f kz=%.2f',Re,kx,kz);
title(str)
filename = sprintf("sigma_with_c/Re=%d_kx=%.2f_kz=%.2f_sigma_with_c.png",Re,kx,kz);
saveas(gcf,filename,'png')