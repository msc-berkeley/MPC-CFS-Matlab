function [ k, K, sigs ] = stochasticLQR( L2, L1, A, B, x1 )
% Compute the LQR step with sigma
nx=length(x1);
nu=size(B,2);
N=size(A,3)+1;
Vxx=L2(1:nx,1:nx,N);
Vx=L1(1:nx,N);
k=zeros(nu,N-1);
K=zeros(nu,nx,N-1);
sigs=zeros(nu,nu,N-1);
for t=N-1:-1:1
    Qxx=L2(1:nx,1:nx,t)+A(:,:,t)'*Vxx*A(:,:,t);
    Quu=L2(nx+1:end,nx+1:end,t)+B(:,:,t)'*Vxx*B(:,:,t);
    Qux=L2(nx+1:end,1:nx,t)+B(:,:,t)'*Vxx*A(:,:,t);
    Qx=L1(1:nx,t)+A(:,:,t)'*Vx;
    Qu=L1(nx+1:end,t)+B(:,:,t)'*Vx;
    sig=inv(Quu);
    sigs(:,:,t)=sig;
    Vx=Qx-Qux'*sig*Qu;
    Vxx=Qxx-Qux'*sig*Qux;
    k(:,t)=-sig*Qu;
    K(:,:,t)=-sig*Qux;
end

end

