function [] = plotResolventModes(y, u, s, v, N, Nsvd, field)

figure
semilogy(s,'o');
ylabel('sigma'); 
xlabel('index');
title('Singular values (primitive)');
grid on
grid minor
% fname = sprintf('%d-%d-%d-%d-singular_values.png',Re,kx,kz,omega);
% saveas(gcf,fname);

% singular response modes
% field = 2; % u = 0, v = 1, w = 2, p = 3

figure
for i = 1:Nsvd
    subplot(Nsvd/2,2,i)
    plot(y, real(u(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
    ylabel('v');
    xlabel('y');
end

sgtitle('Response modes (primitive)')
% fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% saveas(gcf,fname);

% singular forcing modes
% field = 2; % u = 0, v = 1, w = 2, p = 3

figure
for i = 1:Nsvd
    subplot(Nsvd/2,2,i)
    plot(y, real(v(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
    ylabel('v');
    xlabel('y');
end

sgtitle('Forcing modes')
% fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% saveas(gcf,fname);
