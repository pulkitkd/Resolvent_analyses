function [response, singVal] = getSingularModes()
%singular modes for a given set of triplets

% inputs: Re, triplets, N, nsvd

Re = 44000;
triplets = [[6, 1, 0.66]; [6, -1, 0.66]; [1, 1, 0.66]; [1, -1, 0.66]];
N = 16;
nsvd = 4;
nt = length(triplets);
singVal = zeros(nsvd, nt);
response = zeros(4*N, nsvd*nt);

for i = 1:nt
    k = triplets(i,1);
    n = triplets(i,2);
    c = triplets(i,3);
    [r,su,ss,sv,U0,yP,UP,dU0,dr,H] = resolventSVD(Re, k, n, c, N, nsvd);
    singVal(:, i) = ss;
    for j = 1:nsvd
        response(:, (i-1)*nsvd + j) = su(:,j);
    end
end

end

