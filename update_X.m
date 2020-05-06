function [X,t] = update_X(T, Omega, alpha,gamma, factor, X)

CC=tic;
dim = size(T);
N = ndims(T);
for i = 1:N
       tmp(i)=alpha(i)/gamma;
end


dim = size(T);
M = cell(ndims(T), 1);
tau = tmp;
normT = norm(T(:));
    tau = tau*factor;
    Msum = 0;
    for i = 1:ndims(T)
        M{i} = Fold(Pro2TraceNorm(Unfold(X, dim, i), tau(i)), dim, i);
        Msum = Msum +  M{i};
    end
    Xlast = X;
    X = Msum / N;
    X(Omega) = T(Omega);
        t=toc(CC);
end
