%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ u, x ] = ILQR_tur( x0, u0, xref, Ts )
N=length(u0)+1;
x=zeros(5,N);
x(:,1)=x0;
u=u0;
du=zeros(2,N-1);
up=zeros(2,N-1);
dx=zeros(5,N);
P=zeros(5,5,N);
b=zeros(5,N);
c=zeros(1,N);
d=zeros(2,N-1);
e=zeros(1,N-1);
A=zeros(5,5,N-1);
B=zeros(5,2,N-1);
invB=zeros(2,2,N-1);
bss=zeros(2,N-1);
c1=1000;
Q=[c1 0 0 0 0;0 c1 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
Qf=Q;%[10 0 0 0;0 10 0 0;0 0 0 0;0 0 0 0];
R=[10 0;0 20];
P(:,:,N)=Qf;
cost=1e10;

iter=10;
thres=0.01;
J=zeros(iter,1);
lineSearchScale=0.5;
%tic;
for j=1:iter
%% Forward Simulation
    for i=2:N
        thk=x(3,i-1);vk=x(4,i-1);thdk=x(5,i-1);
        ak=u(1,i-1);thak=u(2,i-1);
        thm=thk+0.5*Ts*thdk;
        
        x(1,i) = x(1,i-1) + vk*Ts*cos(thm);
        x(2,i) = x(2,i-1) + vk*Ts*sin(thm);        
        x(3,i) = x(3,i-1) + Ts*thdk;
        x(4,i) = x(4,i-1) + Ts*ak;
        x(5,i) = x(5,i-1) + Ts*thak;
        
%         vk=x(3,i-1);thk=x(4,i-1);ak=u(1,i-1);curk=u(2,i-1);
%         lk=vk*Ts+1/2*Ts^2*ak;thm=thk+curk*lk;
%         if curk~=0
%             x(1,i)=x(1,i-1)+(sin(thm)-sin(thk))/curk;
%             x(2,i)=x(2,i-1)+(cos(thk)-cos(thm))/curk;
%         else
%             x(1,i)=x(1,i-1)+lk*cos(thk);
%             x(2,i)=x(2,i-1)+lk*sin(thk);
%         end
%         x(3,i)=vk+Ts*ak;
%         x(4,i)=thk+curk*lk;
    end
%     plot(x(1,:),x(2,:),'.');
%     pause(1);
    
%% LQR

    xd=xref-x;
    b(:,N)=-Qf*xd(:,N);
    c(N)=1/2*xd(:,N)'*Qf*xd(:,N);
    for k=N:-1:2
        ud=-u;
        thk=x(3,i-1);vk=x(4,i-1);thdk=x(5,i-1);
        ak=u(1,i-1);thak=u(2,i-1);
        thm=thk+0.5*Ts*thdk;
        
%         vk=x(3,k-1);thk=x(4,k-1);ak=u(1,k-1);curk=u(2,k-1);
%         lk=vk*Ts+1/2*Ts^2*ak;thm=thk+curk*lk;

        Ak = [ 1  0  -vk*Ts*sin(thm) Ts*cos(thm) -0.5*Ts^2*vk*sin(thm);
               0  1   vk*Ts*cos(thm) Ts*sin(thm)  0.5*Ts^2*vk*cos(thm);
               0  0         1            0                Ts          ;
               0  0         0            1                 0          ;
               0  0         0            0                 1          ];
        
        Bk = [ 0  0 ;
               0  0 ;
               0  0 ;
               Ts 0 ;
               0  Ts];
           
%         if curk~=0
%             Ak=[1 0 cos(thm)*Ts (cos(thm)-cos(thk))/curk;
%                       0 1 sin(thm)*Ts (sin(thm)-sin(thk))/curk;
%                         0 0 1 0;0 0 curk*Ts 1];
%             Bk=[1/2*Ts^2*cos(thm) (cos(thm)*lk*curk-sin(thm)+sin(thk))/curk^2;
%                 1/2*Ts^2*sin(thm) (sin(thm)*lk*curk+cos(thm)-cos(thk))/curk^2;
%                 Ts 0;
%                 1/2*curk*Ts^2 lk];
%         else
%             Ak=[1 0 Ts*cos(thk) -lk*sin(thk);
%                         0 1 Ts*sin(thk) lk*cos(thk);
%                         0 0 1 0;0 0 curk*Ts 1];
%             Bk=[cos(thk)*Ts^2/2 0;sin(thk)*Ts^2/2 0;Ts 0;1/2*curk*Ts^2 lk];
%         end
        A(:,:,k-1)=Ak;
        B(:,:,k-1)=Bk;
        Pk=P(:,:,k);
        d(:,k-1)=-R*ud(:,k-1);
        e(k-1)=1/2*ud(:,k-1)'*R*ud(:,k-1);
        
        VB=Bk'*Pk*Bk+R;
        bs=Bk'*b(:,k)+d(:,k-1); %%
        inVB=inv(VB);
        P(:,:,k-1)=-Ak'*Pk*Bk*inVB*Bk'*Pk*Ak+Ak'*Pk*Ak+Q;
        b(:,k-1)=-Ak'*Pk*Bk*inVB*bs+Ak'*b(:,k)-Q*xd(:,k-1);
        c(k-1)=c(k)+e(k-1)+1/2*xd(:,k-1)'*Q*xd(:,k-1)-1/2*bs'*inVB*bs;
        invB(:,:,k-1)=inVB; 
        bss(:,k-1)=bs;
    end
    %% Line Search
    flag=0;
    step=1;
    while flag==0
        for k=1:N-1
            du(:,k)=-invB(:,:,k)*(B(:,:,k)'*P(:,:,k+1)*A(:,:,k)*dx(:,k)+step*bss(:,k));
            dx(:,k+1)=A(:,:,k)*dx(:,k)+B(:,:,k)*du(:,k);
        end
        up=u+du;
    %% Forward Simulation and Cost Evaluation
        costs=0;
        for i=2:N
            thk=x(3,i-1);vk=x(4,i-1);thdk=x(5,i-1);
            ak=u(1,i-1);thak=u(2,i-1);
            thm=thk+0.5*Ts*thdk;
%             vk=x(3,i-1);thk=x(4,i-1);ak=up(1,i-1);curk=up(2,i-1);
%             lk=vk*Ts+1/2*Ts^2*ak;thm=thk+curk*lk;

            x(1,i) = x(1,i-1) + vk*Ts*cos(thm);
            x(2,i) = x(2,i-1) + vk*Ts*sin(thm);        
            x(3,i) = x(3,i-1) + Ts*thdk;
            x(4,i) = x(4,i-1) + Ts*ak;
            x(5,i) = x(5,i-1) + Ts*thak;
%             if curk~=0
%                 x(1,i)=x(1,i-1)+(sin(thm)-sin(thk))/curk;
%                 x(2,i)=x(2,i-1)+(cos(thk)-cos(thm))/curk;
%             else
%                 x(1,i)=x(1,i-1)+lk*cos(thk);
%                 x(2,i)=x(2,i-1)+lk*sin(thk);
%             end
%             x(3,i)=vk+Ts*ak;
%             x(4,i)=thk+curk*lk;
            % The costs
            if i<N
                costs=costs+1/2*(x(:,i)-xref(:,i))'*Q*(x(:,i)-xref(:,i))+1/2*up(:,i)'*R*up(:,i);
            else
                costs=costs+1/2*(x(:,i)-xref(:,i))'*Qf*(x(:,i)-xref(:,i));
            end         
        end
        % jump out if the cost decrease
        if costs<=cost||norm(du,2)/norm(u,2)<thres%||step<0.3
            cost=costs;
            flag=1;
        end
        step=step*lineSearchScale;
    end
    u=up;
    J(j)=cost;
    
    if j>1 
        if norm(du,2)/norm(u,2)<thres%||abs(J(j)-J(j-1))/J(j)<thres
            break;
        end
    end

end
%toc;
%plot(x(1,:),x(2,:),'.');

end

