function [ x,cost ] = forwardSim_car( L2,L1,l,x1,u,var,xref )
Ts=var.dt;
N=var.N;
nx = var.nx;
x=zeros(nx,N);
x(:,1)=x1;
cost=0;



for i=2:N
%     thk=x(4,i);curk=u(2,i);lk=x(3,i)*dt+0.5*dt^2*u(1,i);th=thk+curk*lk;
%     if curk>=0.001
%         x(1,i+1)=x(1,i)+(sin(th)-sin(thk))/curk;
%         x(2,i+1)=x(2,i)+(cos(thk)-cos(th))/curk;
%         x(3,i+1)=x(3,i)+dt*u(1,i);
%         x(4,i+1)=th;
%     else
%         x(1,i+1)=x(1,i)+lk*cos(thk);
%         x(2,i+1)=x(2,i)+lk*sin(thk);
%         x(3,i+1)=x(3,i)+dt*u(1,i);
%         x(4,i+1)=th;
%     end
        thk=x(3,i-1);vk=x(4,i-1);thdk=x(5,i-1);
        ak=u(1,i-1);thak=u(2,i-1);
        thm=thk+0.5*Ts*thdk;
        
        x(1,i) = x(1,i-1) + vk*Ts*cos(thm);
        x(2,i) = x(2,i-1) + vk*Ts*sin(thm);        
        x(3,i) = x(3,i-1) + Ts*thdk;
        x(4,i) = x(4,i-1) + Ts*ak;
        x(5,i) = x(5,i-1) + Ts*thak;
        
%         x_f = [x_f;x(:,i)];
%         u_f = [u_f; u(:,i-1)];
        dz = [x(:,i) - xref(:,i);-u(:,i-1)];
        cost = cost + 0.5*dz'*L2(:,:,i)*dz;
        
end
 
% [ a,d ] = findpt( x,ref );
% [ d, a] = derivative( x(1:2,:),ref );
% % for cost of velocity
% wv=var.wv;
% vr=var.vr;
% as=[0 1;-1 0]*a;
% for i=1:N-1
%     cost=cost+0.5*[x(:,i);u(:,i)]'*L2(:,:,i)*[x(:,i);u(:,i)]+[x(:,i);u(:,i)]'*L1(:,i)+l(i);
%     cost=cost+wcur*x(3,i)^2*(sqrt(u(2,i)^2+alpha^2)-alpha); % Cost for curvature
%     if ~var.pure
%         % cost of velocity
%         A=as(:,i)*as(:,i)';
%         B=as(:,i)*vr;
%         C=0.5*vr^2;
%         th=[cos(x(4,i));sin(x(4,i))];
%         cost=cost+wv*(0.5*x(3,i)^2*th'*A*th-x(3,i)*th'*B+C);
%     else
%         % cost of pure velocity tracking
%         cost=cost+0.5*wv*(x(3,i)-vr)^2;
%     end
% end
% cost=cost+0.5*x(:,N)'*L2(1:nx,1:nx,N)*x(:,N)+x(:,N)'*L1(1:nx,N)+l(N);
% cost=cost+0.5*var.wref*sum(d.^2);
% if ~var.pure
%     % cost of velocity
%     A=as(:,N)*as(:,N)';
%     B=as(:,N)*vr;
%     C=0.5*vr^2;
%     th=[cos(x(4,N));sin(x(4,N))];
%     cost=cost+wv*(0.5*x(3,N)^2*th'*A*th-x(3,N)*th'*B+C);
% else
%     % cost of pure velocity tracking
%     cost=cost+0.5*wv*(x(3,N)-vr)^2;
% end

end

