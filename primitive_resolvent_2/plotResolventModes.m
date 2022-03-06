function [] = plotResolventModes(y, u, s, v, N, Nsvd, field, yc, Re, kx, kz, omega)

figure(1)
% subplot(1,2,1)
semilogy(s,'o','linewidth',2);
ylabel('sigma'); 
xlabel('index');
title('Singular values (primitive)');
grid on
grid minor
fname = sprintf('%d-%d-%d-%d-%d-singular_values.png',Re,kx,kz,omega,N);
saveas(gcf,fname);
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
        plot(y, imag(u(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
        xline(yc,'r')
        ylabel('v');
        xlabel('y');
    end
else
% If field ~= 0,1,2 plot all fields
    for j = 0:2
        for i = 1:Nsvd
            subplot(Nsvd/2,2,i)
            plot(y, imag(u(j*N+1:(j+1)*N,i)), 'LineWidth', 2)
            xline(yc,'r')
            ylabel('v');
            xlabel('y');
            hold on
        end
    end
end
sgtitle('Response modes')
% hold off
fname = sprintf('%d-%d-%d-%d-%d-response_modes.png',Re,kx,kz,omega,N);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,fname);

% singular forcing modes
% field = 2; % u = 0, v = 1, w = 2, p = 3

figure(3)
if (field == 0 || field == 1 || field == 2)
    for i = 1:Nsvd
        subplot(Nsvd/2,2,i)
        plot(y, imag(v(field*N+1:(field+1)*N,i)), 'LineWidth', 2)
        xline(yc,'r')
        ylabel('uvw');
        xlabel('y');
    end
else
% If field ~= 0,1,2 plot all fields
    for j = 0:2
        for i = 1:Nsvd
            subplot(Nsvd/2,2,i)
            plot(y, imag(v(j*N+1:(j+1)*N,i)), 'LineWidth', 2)
            xline(yc,'r')
            ylabel('u,v,w');
            xlabel('y');
            hold on
        end
    end
end
sgtitle('Forcing modes')
hold off
fname = sprintf('%d-%d-%d-%d-%d-forcing_modes.png',Re,kx,kz,omega,N);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,fname);
