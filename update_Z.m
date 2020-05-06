function [X,t, U, Theta, V, numiter ] = update_Z(size_tens, r, Known, data, opts,X,beta1 )
BB= tic;
s1 = size_tens(1);
s2 = size_tens(2);
s3 = size_tens(3);
m = s1*s2;
n = s3;

 
 if nargin < 5
    epsilon = 1e-5; % convergence threshold
else
    epsilon = opts.epsilon;
end
X=reshape(X,m,n);
beta = opts.beta;
X(Known)=data;
[indm,indn,data1]=find(X);
data( data1 == 0 )= eps;
res = sparse(indm, indn, data1, m, n);
[indm, indn, data1] = find(res);
 
Msum = [];

verbosity = 1;
printsyb = ['-', 'X', '|'];

i = 1;
W = 0;
 
U = []; V = []; Theta = [];
nnorm = norm(data1, 'fro');

gamma = 1;

 x3_old = ones(n,1);
 x3_old = x3_old/norm(x3_old);
 
while (i <= r) && (gamma > epsilon )
        i = i + 1;
    % 1. find the top singular pair of the residual and update the gresnorm
    resvec = beta1*(data1 - W);      % compact vector representation of the negative gradient
    sparse_update(res, resvec); % sparse update the res using resvec
     [u,v] = get_init2(-res,s1,s2,x3_old);
    x3_old = v;
    % 2. update the weight Theta, the pursuit basis is uv', its weight is s.
    S = -beta*sparse_inp(u', v', indm, indn)';       % the compact vector representation of u*v'
    d = S - W;
 
    
    tmp =  (dot(resvec,d));
    gamma =  abs(tmp)*1.2/norm(d)^2;              % faster than the recursively updating rule...    
 
    U = U*(1-gamma);
    U = [U -gamma*beta*u];
    V = [V v];
     
    W = W + gamma*d;   
end
% V = diag(Theta)*V;
X = U*V';
% if opts.no_noise == 1
%     X(Known) = data;
% end
numiter = i;
 t=toc(BB);

 


function [vu,w,e] = get_init2(mx,s1,s2,x3_old)

%     mx = tens2mat(X,1:2,3);
%     mmx = mx'*mx;
%     [w,e] = eigs(mmx,1,'lm');

    [w,e]=power_method(mx,x3_old,2);

    v = mx*w/sqrt(e);
    
    
%     mv = vec2tens(v,[size(X,2) size(X,3)],1:2);
    improved_times = 1;
    for ii = 1: improved_times
        mv = reshape(v,s1, s2);
        mmv = mv*mv';

        [u,e1] = eigs(mmv,1,'lm');
        v = mv'*u/sqrt(e1);

        vu = kron(v,u);

        w =  vu'*mx;          % kron is opposed to the direction of tens2mat
        w = w'/norm(w);
    end
