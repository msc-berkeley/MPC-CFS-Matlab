function [ x,u,xbar,ubar,k1,K,sigs, A, B ] = MaxEntILQR( L2, L1, l, x0, u, var, xref )
dt=var.dt;
[x,cost]=forwardSim_car( L2,L1,l,x0,u,var,xref );
N = var.N;
maxstep=100;
J=zeros(maxstep,1);
J(1)=cost;
Thres=var.Thres;
lineSearchThres=var.lineSearchThres;
% plot(x(1,:),x(2,:),'.');hold on;
% plot(ref(1,:),ref(2,:));
% axis([-5 30 -5 50]);
% axis equal;hold off;
% pause(0.5);
for i=1:maxstep
    i;
    xbar=x;
    ubar=u;
    stopsearch=0;
    [ A,B ] = linearizeAB_car( x(:,1:end-1),u,var ); 
    L1 = modifyCost( L2, L1,x,u,xref ); % Transform cost from traj to delta traj 
    %[dL2,dL1] = nonlinearCost(x,u,var,ref); % Add the non-quadratic cost function
    
    [ k, K, sigs ] = stochasticLQR( L2, L1, A, B, x0 );
    k1=k;
    % Line search
    step=1;
    while 1
        for t=1:N-1
            k1(:,t)=u(:,t)+step*k(:,t)-K(:,:,t)*x(:,t);% total control
        end
        [ x1,u1,cost1 ] = executePolicy( L2,L1,l,x0,k1,K,var,xref );
        if cost1<cost
            cost=cost1;
            J(i+1)=cost;
            x=x1;
            u=u1;
            k=k1;
            break;
        end
        step=step*0.5;
        if step<lineSearchThres
            stopsearch=1;
            x=x1;
            u=u1;
            break;
        end
    end
    if stopsearch==1
        break;
    end
    if (J(i)-J(i+1))/J(i+1)<Thres
        break;
    end
%     plot(x(1,:),x(2,:),'.');hold on;
%     plot(xref(1,:),xref(2,:));
%     axis([-1 8 -2 2]);
%     axis equal;hold off;
%     pause(0.05);
end

end

