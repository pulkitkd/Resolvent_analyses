%program to evaluate interaction coefficients for interaction between two
%Fourier modes - K1 (variable) and K2 (fixed)

clear
clc

itri = 15;
itrj = 15;
intCoeff = zeros(itri,itrj);

Re = 187;
N  = 100;
Nsvd = 10;
[y, ~, ~, ~, U0] = channelMeanVel(Re, N);

[~,DM] = chebdif(N,2); 
D1 = DM(:,:,1); % first derivative - N X N
Z = zeros(N);

[W, ~] = weight_matrix(N);
W = W(1:3*N,1:3*N);
% test=u'*(W'*W)*u;

for i = 1:itri
    for j = 1:itrj
        K1 = [i            ,j            , 1];
        K2 = [0            ,1            , 1];
        K3 = [K1(1) + K2(1),K1(2) + K2(2), 1];
        disp(K1)
        
        kx   = [K1(1), K2(1), K3(1)];
        kz   = [K1(2), K2(2), K3(2)];
        c    = [K1(3), K2(3), K3(3)];
        omega= c.*kx;
        
        %% get the resolvent modes
        [psi1, s1, phi1] = getResolventSVD(kx(1), kz(1), omega(1), Re, N, Nsvd, U0);
        [psi1, ~] = channelSym(psi1, phi1, s1, Nsvd);
        
        [psi2, s2, phi2] = getResolventSVD(kx(2), kz(2), omega(2), Re, N, Nsvd, U0);
        [psi2, ~] = channelSym(psi2, phi2, s2, Nsvd);
        
        [psi3, s3, phi3] = getResolventSVD(kx(3), kz(3), omega(3), Re, N, Nsvd, U0);
        [psi3, phi3] = channelSym(psi3, phi3, s3, Nsvd);
        % plotResolventModes(y, psi2, s2, phi2, N, Nsvd, 0);
        
        %% compute the interaction coeff
        % Ensure that we are picking the right even and odd modes
        % v should be odd, u and w even
        % when we plot all three fields, it is seen that 
        % -first column of u and w is even
        % -first column of v is odd
        % ...conveniently ?
        
        u1hat = iseven(psi1(1:N,1),psi1(1:N,2));
        v1hat = isodd(psi1(N+1:2*N,1),psi1(N+1:2*N,2));
        w1hat = iseven(psi1(2*N+1:3*N,1),psi1(2*N+1:3*N,2));
        
        f3uhat = iseven(phi3(1:N,1),phi3(1:N,2));
        f3vhat = isodd(phi3(N+1:2*N,1),phi3(N+1:2*N,2));
        f3what = iseven(phi3(2*N+1:3*N,1),phi3(2*N+1:3*N,2));
        f3hat = [f3uhat; f3vhat; f3what];
        
        u2hat = iseven(psi2(1:N,1),psi2(1:N,2));
        v2hat = isodd(psi2(N+1:2*N,1),psi2(N+1:2*N,2));
        w2hat = iseven(psi2(2*N+1:3*N,1),psi2(2*N+1:3*N,2));
        % 
        % % Check that u1hat , w1hat are even and v1hat is odd
        % % Plot the modes if they aren't even/odd so the user can inspect visually
        % if (isevenodd(u1hat)~=0 || isevenodd(w1hat)~=0 || isevenodd(v1hat)~=1)
        %     disp('u v w are not appropriately even or odd')
        %     plot(y,real(u1hat))
        %     hold on
        %     plot(y,real(v1hat))
        %     plot(y,real(w1hat))
        %     hold off
        % end
        % 
        % % Check that u2hat , w2hat are even and v2hat is odd
        % if (isevenodd(u2hat)~=0 || isevenodd(w2hat)~=0 || isevenodd(v2hat)~=1)
        %     disp('u v w are not appropriately even or odd')
        %     plot(y,real(u2hat))
        %     hold on
        %     plot(y,real(v2hat))
        %     plot(y,real(w2hat))
        %     hold off
        % end
        
       
        %% compute the interaction coefficient
        
        divpsipsi =        1i*kx(3)*[u1hat.*u2hat; u1hat.*v2hat; u1hat.*w2hat]+...
           [D1,Z,Z; Z,D1,Z; Z,Z,D1]*[v1hat.*u2hat; v1hat.*v2hat; v1hat.*w2hat]+...
                           1i*kz(3)*[w1hat.*u2hat; w1hat.*v2hat; w1hat.*w2hat];

        innerprod = divpsipsi'*(W'*W)*f3hat;
        intCoeff1 = -s1(1)*s2(1)*innerprod;
        
        divpsipsi =          1i*kx(3)*[u1hat.*u2hat; v1hat.*u2hat; w1hat.*u2hat]+...
             [D1,Z,Z; Z,D1,Z; Z,Z,D1]*[u1hat.*v2hat; v1hat.*v2hat; w1hat.*v2hat]+...
                             1i*kz(3)*[u1hat.*w2hat; v1hat.*w2hat; w1hat.*w2hat];
        
        innerprod = divpsipsi'*(W'*W)*f3hat;
        intCoeff2 = -s1(1)*s2(1)*innerprod;
        
        intCoeff(i,j) = (intCoeff1 + intCoeff2);

        if (K2(1)==0 && K2(2)==0)
            intCoeff(i,j)=NaN;
        end
        
        singVals = s3(1);
        responseModes = psi3(:,1);
    end
end

%%
figure(1)
heatmap(real(intCoeff))
colorbar('eastoutside')
xlabel('kz')
ylabel('kx')
title('real part of interaction coefficient')

figure(2)
heatmap(imag(intCoeff))
colorbar('eastoutside')
xlabel('kz')
ylabel('kx')
title('imag part of interaction coefficient')

figure(3)
heatmap(abs(real(intCoeff)./imag(intCoeff)))
colorbar('eastoutside')
xlabel('kz')
ylabel('kx')
title('ratio of real to imag parts')
% figure(2)
% surf(abs(real(intCoeff)))
% colorbar('eastoutside')
%%
% Important for GQL ??!
%%
% hold on
% plot(real(intCoeff),'r-*')
% hold off

%%

% surf(abs(imag(intCoeff))/abs(real(intCoeff)))
% colorbar('eastoutside')
