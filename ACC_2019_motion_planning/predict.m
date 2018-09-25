function [ V0, mus, sigs ] = predict( L2, L1, l, A, B, x1, Xn2, xbar, ubar, k, K, sigsu )
% Calculate the mus and sigs of the predicted trajectory; Calculate the
% cost-to-go V0; Draw the ellipse
% [ k, K, sigsu, V0 ] = MaxEntLQR( L2, L1, l, A, B, x1 );
N=size(A,3)+1;
nx=length(x1);
nu=size(B,2);
x=zeros(nx,N);
x(:,1)=x1;
sigs=zeros(4,4,N);
sigs(:,:,1)=0.01*eye(4);
V=[1 0 0 0;0 1 0 0];
for t=1:N-1
    u=k(:,t)+K(:,:,t)*x(:,t);
    x(:,t+1)=A(:,:,t)*x(:,t)+B(:,:,t)*u+xbar(:,t+1)-A(:,:,t)*xbar(:,t)-B(:,:,t)*ubar(:,t);
    sigs(:,:,t+1)=(A(:,:,t)+B(:,:,t)*K(:,:,t))*sigs(:,:,t)*(A(:,:,t)+B(:,:,t)*K(:,:,t))'+B(:,:,t)*sigsu(:,:,t)*B(:,:,t)';
    ellipse(x(1:2,t),Xn2*V*sigs(:,:,t)*V');hold on;
end
% axis equal;
% axis([-5 30 -2 30]);
mus=x(1:2,:);
end

