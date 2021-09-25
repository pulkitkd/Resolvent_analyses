function Q = generateWeightMatrix(N,k2)

% Weight matrix for energy inner product


[~,DM] = chebdif(N,2);
D1T = DM(:,:,1);
[~,dy] = clencurt(N-1);
IWT = diag(dy);
IW = IWT(2:end-1, 2:end-1);
QvT = (D1T'*IWT*D1T/k2 + IWT);
Qv = QvT(2:end-1,2:end-1);
Qeta = IW/k2;
Z = zeros(N-2);
Q = [Qv, Z; Z, Qeta];
