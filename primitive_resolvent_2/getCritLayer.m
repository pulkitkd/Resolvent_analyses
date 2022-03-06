function [yc] = getCritLayer(y,U0,c)
% Critical layer height
Y = -1:0.01:1;
fU0 = griddedInterpolant(flip(y),flip(U0));
yc = 0;

for i=2:length(Y)-1
    if (fU0(Y(i))-c)*(fU0(Y(i+1))-c) < 0
        yc = Y(i);
    end
end

end