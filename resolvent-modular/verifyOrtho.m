[y,w] = clencurt(N-1);
W = blkdiag(diag(w),diag(w),diag(w),diag(w));
I = eye(length(y));
Z = zeros(length(y));
M = [I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z];
W = W*M;
dot(usu(:,1), usu(:,2))