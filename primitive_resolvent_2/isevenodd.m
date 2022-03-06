function [a] = isevenodd(u)

tol = 0.01*max(abs(u));
u = real(u);
    if tol < 1e-15
        disp('Possible zero mode')
        a = -1;
    elseif norm(u - flip(u)) < tol
        %function is even
        a = 0;
    elseif norm(u + flip(u)) < tol
        %function is odd
        a = 1;
    else
        %function is neither even nor odd
        disp('Error! Function is neither even nor odd.')
        a = 2;
        disp(tol)
        disp(norm(u - flip(u)))
        disp(norm(u + flip(u)))
    end
end