function [a] = isevenodd(u)

tol = 1e-5;
u = real(u);
    if norm(u - flip(u)) < tol
        %function is even
        a = 0;
    elseif norm(u + flip(u)) < tol
        %function is odd
        a = 1;
    else
        %function is neither even nor odd
        a = 2;
    end
end