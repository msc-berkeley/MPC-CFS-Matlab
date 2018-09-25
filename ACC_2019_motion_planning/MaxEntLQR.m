function [ k, K, sigs, V0 ] = MaxEntLQR( L2, L1, l, A, B, x1 )
% Solve the maximum entropy LQR problem with cost of
% [x;u]'L2[x;u]+[x;u]'L1+l; and linear system x(t+1)=A(t)x(t)+B(t)u(t); and
% initial trajectory x1. The output is the policy N(k(t)+K(t)x(t);sig(t)) and
% the value of initial state V0.
nx=length(x1);
nu=size(B,2);
N=size(A,3)+1;
P=L2(1:nx,1:nx,N);
q=L1(1:nx,N);
b=l(N);
k=zeros(nu,N-1);
K=zeros(nu,nx,N-1);
sigs=zeros(nu,nu,N-1);
for t=N-1:-1:1
    % The intermedium variable
    sig=inv(L2(nx+1:end,nx+1:end,t)+B(:,:,t)'*P*B(:,:,t));
    sigs(:,:,t)=sig;
    D=L2(1:nx,1:nx,t)+A(:,:,t)'*P*A(:,:,t);
    E=L2(nx+1:end,1:nx,t)+B(:,:,t)'*P*A(:,:,t);
    F=L1(1:nx,t)+A(:,:,t)'*q;
    G=L1(nx+1:end,t)+B(:,:,t)'*q;
    % Update V
    P=D-E'*sig*E;
    q=F-E'*sig*G;
    b=l(t)+b-1/2*G'*sig*G-2*pi*sqrt(det(sig));
    % The policy
    k(:,t)=-sig*G;
    K(:,:,t)=-sig*E;
end
V0=1/2*x1'*P*x1+x1'*q+b;

