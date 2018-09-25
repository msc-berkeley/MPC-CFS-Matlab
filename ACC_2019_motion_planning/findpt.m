function [ a,d ] = findpt( x,ref )
% given a polyline ref and a traj, find the corresponding vector for each
% point: a, and the corrresponding distance: d
N=size(x,2);
seg=ref(:,2:end)-ref(:,1:end-1);
as=[0 -1;1 0]*seg;
as=normc(as);
% dis=abs(as'*x(1:2,:)-repmat(diag(as'*ref(:,1:end-1)),1,N));
% [d,num]=min(dis);
% a=as(:,num);
a=zeros(2,N);
d=zeros(N,1);
for i=1:N
%     disvec=x(1:2,i)*ones(1,N)-ref;
%     sum((x(1:2,i)*ones(1,size(ref,2))-ref).^2)
    [val, num]=min(sum((x(1:2,i)*ones(1,size(ref,2))-ref).^2));
    x0=ref(:,num);xpre=ref(:,num-1);xpos=ref(:,num+1);
    th1=acos(((xpre-x0)/norm(xpre-x0))'*((x(1:2,i)-x0)/norm(x(1:2,i)-x0)));
    th2=acos(((xpos-x0)/norm(xpos-x0))'*((x(1:2,i)-x0)/norm(x(1:2,i)-x0)));
    if th1<th2
        a(:,i)=as(:,num-1);
        d(i)=abs(a(:,i)'*(x(1:2,i)-ref(:,num-1)));
    else
        a(:,i)=as(:,num);
        d(i)=abs(a(:,i)'*(x(1:2,i)-ref(:,num)));
    end
end

end

