clc;
clear;
filename='TestImages/peppers.bmp';

img = imread(filename);
T= im2double(img);
size_tensor=size(T);
Nway =size_tensor;
ObsRatio=20;                 
Omega = randperm(prod(Nway));
Omega = Omega(1:round((ObsRatio/100)*prod(Nway)));
O = zeros(Nway);
O(Omega) = 1;
y=T.*O;

% y1(:,:,:,ObsRatio)=y;
known=find(y);
data=y(known);
multi1=zeros(Nway);
multi2=zeros(Nway);
r = 20;
opts.no_noise = 1;
opts.epsilon = 1e-5;
opts.beta = 6500;
alpha = [1, 1, 1e-3];
alpha = alpha / sum(alpha);
maxIter = 50;
epsilon = 1e-5;
beta = 0.7;
stop_iter=1;
X=y;
Lambda=zeros(size_tensor);
Xlast=zeros(size_tensor);
kkk=0;
for i=1:maxIter
     X=X+Lambda/beta;
     known=find(X);
      data=X(known);
     AA=tic;
     [Z,t] =update_Z(size_tensor, r,known, data, opts,X,beta);
     CP(i)=t;
      Z=reshape(Z,size_tensor);
      lambda = 0.99;
      Z=Z-Lambda/beta;
      [X,t] = update_X( T,Omega,alpha, beta,lambda,Z );
      TT(i)=t;
  %    Lambda=zeros(size_tensor);
      Lambda=Lambda+beta*(X-Z);
      ADMM(i)=toc(AA);
      relative_error(i)=norm(X(:)-Xlast(:),'fro')/norm(T(:),'fro');
     beta=beta*1.05; 
if(relative_error(i)<1e-5)
    break;
end
%  kkk=relative_error(i);
 Xlast = X;
end
figure;
subplot(2,2,1),imshow(T),title('orginal image');
subplot(2,2,2),imshow(y),title('missing image');
subplot(2,2,3),imshow(X),title('recovered image');
subplot(2,2,4),
i=2:1:length(relative_error);
plot(i,relative_error(i)),title('check convergence'),xlabel('iteration'),ylabel('relative error');
%   lambda=lambda+5;

     