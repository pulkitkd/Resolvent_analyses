function [c] = iseven(a,b)
% returns an even function out of a and b

if isevenodd(a) == -1
    disp('Found possible zero mode, setting it zero')
    c = zeros(length(a),1);
elseif isevenodd(a) == 0
    c = a;
elseif isevenodd(b) == 0
    c = b;
else
    disp("None of the functions are even. Error is")
    disp(norm(a - flip(a)))
    disp(norm(b - flip(b)))
    c = zeros(length(a),1);
end

end