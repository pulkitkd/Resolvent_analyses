function [W, iW] = weight_matrix(N)

[~,w] = clencurt(N-1);

w = sqrtm(diag(w));

Z = zeros(N);
I = eye(N);
iw = I/w;

W  = [w  Z Z Z; Z w  Z Z; Z Z w  Z; Z Z Z Z];
iW = [iw Z Z Z; Z iw Z Z; Z Z iw Z; Z Z Z Z];

end