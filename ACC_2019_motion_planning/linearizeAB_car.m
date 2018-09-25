function [ A,B ] = linearizeAB_car( x,u,var )
Ts=var.dt;
N=var.N;
nx=var.nx;
nu=var.nu;
A=zeros(nx,nx,N-1);
B=zeros(nx,nu,N-1);

for i=2:N
   
     thk=x(3,i-1);vk=x(4,i-1);thdk=x(5,i-1);
     ak=u(1,i-1);thak=u(2,i-1);
     thm=thk+0.5*Ts*thdk;
       

     A(:,:,i-1) = [ 1  0  -vk*Ts*sin(thm) Ts*cos(thm) -0.5*Ts^2*vk*sin(thm);
                   0  1   vk*Ts*cos(thm) Ts*sin(thm)  0.5*Ts^2*vk*cos(thm);
                   0  0         1            0                Ts          ;
                   0  0         0            1                 0          ;
                   0  0         0            0                 1          ];
        
     B(:,:,i-1) = [ 0  0 ;
                   0  0 ;
                   0  0 ;
                   Ts 0 ;
                   0  Ts];

    
%     thk=x(4,i); curk=u(2,i);lk=x(3,i)*dt+0.5*u(1,i)*dt^2;th=thk+curk*lk;
%     if curk>=0.001
%         A(:,:,i)=[1 0 cos(th)*dt (cos(th)-cos(thk))/curk;
%             0 1 sin(th)*dt (sin(th)-sin(thk))/curk;
%             0 0 1 0; 0 0 curk*dt 1];
%         B(:,:,i)=[0.5*dt^2*cos(th) (cos(th)*curk*lk-sin(th)+sin(thk))/curk^2;
%             0.5*dt^2*sin(th) (sin(th)*curk*lk+cos(th)-cos(thk))/curk^2;
%             dt 0;0.5*curk*dt^2 lk];
%     else
%         A(:,:,i)=[1 0 cos(thk)*dt -lk*sin(thk);
%             0 1 sin(thk)*dt lk*cos(thk);
%             0 0 1 0;0 0 0 1];
%         B(:,:,i)=[0.5*dt^2*cos(thk) 0;0.5*dt^2*sin(thk) 0;dt 0;0 lk];
%     end
end
end

