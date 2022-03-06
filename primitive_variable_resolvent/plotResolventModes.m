function [] = plotResolventModes(y, u, s, v, N, Nsvd, field)

figure(1)
% subplot(1,2,1)
semilogy(s,'o','linewidth',2);
ylabel('sigma'); 
xlabel('index');
title('Singular values (primitive)');
grid on
grid minor
% fname = sprintf('%d-%d-%d-%d-singular_values.png',Re,kx,kz,omega);
% saveas(gcf,fname);
% subplot(1,2,2)
% semilogy(s./amplitude','o');
% ylabel('sigma/amp'); 
% xlabel('index');
% title('Singular values (scaled by amp of psi modes)');
% grid on
% grid minor




% singular response modes
% field = 0 for u, 1 for v, 2 for w
figure(2)
if (field == 0 || field == 1 || field == 2)
    for i = 1:Nsvd
        subplot(Nsvd/2,2,i)
        plot(y, real(u(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
        ylabel('v');
        xlabel('y');
    end
else
% If field ~= 0,1,2 plot all fields
    for field = 0:2
        for i = 1:Nsvd
            subplot(Nsvd/2,2,i)
            plot(y, real(u(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
            ylabel('v');
            xlabel('y');
        end
    end
end
% sgtitle('Response modes')

% fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% saveas(gcf,fname);

% % singular forcing modes
% % field = 2; % u = 0, v = 1, w = 2, p = 3
% 
% figure
% for i = 1:Nsvd
%     subplot(Nsvd/2,2,i)
%     plot(y, real(v(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
%     ylabel('v');
%     xlabel('y');
% end
% 
% sgtitle('Forcing modes')
% % fname = sprintf('%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% % saveas(gcf,fname);
