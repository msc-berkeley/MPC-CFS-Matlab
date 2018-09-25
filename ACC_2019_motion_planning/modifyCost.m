function [ L1 ] = modifyCost( L2, L1,x,u,xref )
nx=size(x,1);
N=size(x,2);

for i=1:N-1
    L1(:,i)=L2(:,:,i)*[x(:,i);u(:,i)]-L2(:,:,i)*[xref(:,i); 0; 0];
end
L1(1:nx,N)=L2(1:nx,1:nx,N)*x(:,N)-L2(1:nx,1:nx,N)*xref(:,i); 

end

