function [psi_sym, phi_sym] = channelSym(psi, phi, sig, nsvd)

% Enforces symmetry of resolvent forcing and response modes about channel
% centerline (y = 0).
% 
% Original version by Rashad Moarref
% Cleaned-up version by Ryan McMullen


N = size(psi, 1)/3;      % number of Chebyshev points
cntr = floor( (N+1)/2 ); % index closest to channel center

psi_sym = psi;
phi_sym = phi;

for i = 1:2:nsvd
                
    [~, ind] = max( psi(1:cntr, i) ); % location of max amplitude u

    % if maximum is at centerline, then use the previous point
    if ind == (N+1)/2
        ind = ind-1;
    end
    
    % if pairs 'unequal', do nothing (these modes are aleady symmertric)
    if ( sig(i)/sig(i+1) > 1.01 ) || ...
       ( abs(psi(cntr, i) - psi(cntr, i+1)) > 1e-3 )     
    
    % otherwise, enforce symmetry
    else

        psi1 = psi(ind, i) - psi(N+1-ind, i);
        psi2 = psi(ind, i+1) - psi(N+1-ind, i+1);

        r1 = 1/sqrt(1 + abs(psi1/psi2)^2);
        r2 = sqrt(1 - r1^2);

        phase1 = unwrap( angle(psi1) );%phase(psi1);
        phase2 = unwrap( angle(psi2) );%phase(psi2);

        theta2 = 0;
        theta1 = theta2 + phase2 - phase1 + pi;

        a1 = r1*exp(1i*theta1);
        a2 = r2*exp(1i*theta2);
        
        % make first mode in pair even, second mode odd about center
        psi_sym(:, i) = a1*psi(:, i) + a2*psi(:, i+1);
        psi_sym(0*N+(1:cntr), i+1) = psi_sym(0*N+(1:cntr), i);
        psi_sym(1*N+(1:cntr), i+1) = psi_sym(1*N+(1:cntr), i);
        psi_sym(2*N+(1:cntr), i+1) = psi_sym(2*N+(1:cntr), i);
        psi_sym(0*N+(cntr+1:N), i+1) = -psi_sym(0*N+(cntr+1:N), i);
        psi_sym(1*N+(cntr+1:N), i+1) = -psi_sym(1*N+(cntr+1:N), i);
        psi_sym(2*N+(cntr+1:N), i+1) = -psi_sym(2*N+(cntr+1:N), i);
        
        % do the same for phi
        phi_sym(:,i) = a1*phi(:, i) + a2*phi(:, i+1);
        phi_sym(0*N+(1:cntr), i+1) = phi_sym(0*N+(1:cntr), i);
        phi_sym(1*N+(1:cntr), i+1) = phi_sym(1*N+(1:cntr), i);
        phi_sym(2*N+(1:cntr), i+1) = phi_sym(2*N+(1:cntr), i);
        phi_sym(0*N+(cntr+1:N), i+1) = -phi_sym(0*N+(cntr+1:N), i);
        phi_sym(1*N+(cntr+1:N), i+1) = -phi_sym(1*N+(cntr+1:N), i);
        phi_sym(2*N+(cntr+1:N), i+1) = -phi_sym(2*N+(cntr+1:N), i);


    end

    % set phase to zero at peak of (rescaled) u
    val1 = max( psi_sym(1:cntr, i) );
    angle1 = exp(-1i*unwrap( angle(val1) ));%exp(-1i*phase(val1));

    psi_sym(:, i) = psi_sym(:, i)*angle1;
    phi_sym(:, i) = phi_sym(:, i)*angle1;
    
    % do the same for the second mode
    val2 = max( psi_sym(1:cntr, i+1) );
    angle2 = exp(-1i*unwrap( angle(val2) ));%exp(-1i*phase(val2));

    psi_sym(:, i+1) = psi_sym(:, i+1)*angle2;
    phi_sym(:, i+1) = phi_sym(:, i+1)*angle2;

end