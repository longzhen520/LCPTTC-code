% requre: A is m*n, m>n.
function [vec,value]=power_method(A,start,max_iter)
%
%Power method for computing eigenvalues
%
dd=1;
x=start;
x = x/norm(x);
n=10;
iter = 1;
while iter <= max_iter
    iter = iter + 1;
y=A*x;
y = A'*y;
% dd=abs(norm(x)-n);
n=norm(y);
x=y/n;
end
vec=x;
value=n;
