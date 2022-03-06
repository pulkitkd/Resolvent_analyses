function [c] = isodd(a,b)
% returns an odd function out of a and b

if isevenodd(a) == -1
    disp('found possible zero mode, setting it zero')
    c = zeros(length(a),1);
elseif isevenodd(a) == 1
    c = a;
elseif isevenodd(b) == 1
    c = b;
else
    disp("None of the functions are odd")
    c = zeros(length(a),1);
end

end