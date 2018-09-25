% The trajectory include x0 and xpre
% Soft constraint for nstep s and side obs
function [ traj ] = Planner_SoftSide( v0, obs, refpath, dt, x0UTM) % refpath at local coord, obs pos and vel at local position
%% The initials of host vehicle
x0=[0;0]; % initial position
theta=0;
v0vec=[v0;0]; % initial velocity vector
x0pre=[-v0*dt;0];
if v0<0.1
    x0pre=[-0.1*dt;0];
end
refpath=[x0pre refpath];
nstep=size(refpath,2); % the preview horizon
L=2.5;
steerMax=30;
ubar=tan(steerMax/180*pi)/L;
SCFS=0;
%% Initials for the algorithm
dim = 2; % state in 2D position
convergeThres=0.01; % Threshold for convergence check
kmax=5; % Maximum iteration
%% The size of host vehicle
wid=1.1; % half the width
len=2.5; % half the length
dj=sqrt(wid^2+wid^2);
margin=dj;
%% The surrounding vehicles
nobs=size(obs,2);
nobj=nobs;% number of surrounding vehicles we will consider
%% The cost function
% The velocity and acceleration computation matrix
Vdiff = -eye(nstep*dim)+diag(ones(1,(nstep-1)*dim),dim);
Adiff = -Vdiff-diag(ones(1,(nstep-1)*dim),dim)+diag(ones(1,(nstep-2)*dim),dim*2);
Adiff=Adiff(1:(nstep-2)*dim,:)/dt^2;
weight=1;
Adiff=blkdiag([1 0;0 weight],eye(2*(nstep-3)))*Adiff;
% The cost matrix for acceleration
Qa = Adiff'*Adiff;
% The weight for distance and acceleration, curve and soft margin
c = [15,10,0,1000];
% The cost for acceleration
Qabs = c(2)*Qa;
% The cost for reference tracking
Qref = eye(nstep*2);
Qref = c(1)*Qref;
%% The boundary constraint
% define the initial state and velocity
if SCFS==1
    Aeq = zeros(2*dim,nstep*dim+2*nstep); % u_sin, u_cos
else
    Aeq = zeros(2*dim,nstep*dim);
end
Aeq(1:2,1:2) = eye(dim);
Aeq(3:4,3:4) = eye(dim);
Aeq=[Aeq zeros(4,nstep)];

% %% Figures
% figure(1);hold on;
% set(gcf,'Position',[0,200,1000,250], 'color','w');
% axis equal

%% The Planning
% The boundary constraint
beq = [x0pre;x0];
% The iteration
path = [];
path=[path;x0pre;x0];
for i=1:nstep-2
    path=[path;x0+min(x0-x0pre,[1;0]*dt)*i];
%     path=[path;x0+v0*i*dt];
%     path=[path;path(end-1:end)+[vs(i+1);0]*dt];
end
path=path(:);
%%%%%%%%%%%%%%%%%%%%%%
% path=refpath(:);
% Add the soft
path=[path;zeros(nstep,1)]; %x,s
% path=refpath(:);
P = [0 1;-1 0];
% P = -[1 0; 0 1];
refinput = [zeros(1,nstep) 0.5*ones(1,nstep)]; % u_sin, u_cos

% For UTM transformation
R=[cos(x0UTM(3)) -sin(x0UTM(3));sin(x0UTM(3)) cos(x0UTM(3))];
p=x0UTM(1:2);

for k = 1:kmax
    % The constraint
    Lstack = []; Sstack = [];
    for i=1:nstep
        if i>2
            % surrounding vehicle constraints
            for j=1:nobj
%                 poly = obs{j}.x0*ones(1,4)+obs{j}.poly+obs{j}.v*ones(1,4)*dt*(i-2);
                if ~isempty(obs{j})
                    p0=path((i-1)*dim+1:i*dim);% the center point
                    [L,S,d] = d2side(p0,obs{j},dt,i-2);
    %                 % Now consider circle
    %                 margin=5;
    %                 center=obs{j}.x0+obs{j}.v*dt*(i-2);% the center of the obstacle
    %                 p0=path((i-1)*dim+1:i*dim)';% the center point
    %                 [ L, S] = d2Circle( center, p0, margin );
                    if SCFS==1
                        if obs{j}.type==1 % if is pedestrian, then only consider obs avoidance after pedestrian passes
                            if obs{j}.timeL*R*(obs{j}.x0+obs{j}.v*dt*(i-2))+obs{j}.timeL*p+obs{j}.timeS<0 % can start plan to gothrough even have not pass now
%                             if obs{j}.timeL*R*obs{j}.x0+obs{j}.timeL*p+obs{j}.timeS<0 % must wait until pass, then begin to plan go through
                                Lbound=-obs{j}.boundL*R;
                                Sbound=obs{j}.boundL*p+obs{j}.boundS;
                                Lstack=[Lstack;zeros(1,(i-1)*dim) Lbound zeros(1,(nstep-i)*dim) zeros(1,2*nstep) zeros(1,i-1) -1 zeros(1,nstep-i)];
                                Sstack=[Sstack;Sbound-margin];
                            else
                                Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,2*nstep) zeros(1,i-1) -1 zeros(1,nstep-i)]; % x, u, s
                                Sstack = [Sstack;S-margin];
                            end
                        else
                            Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,2*nstep) zeros(1,i-1) -1 zeros(1,nstep-i)];
                            Sstack = [Sstack;S-margin];
                        end
                    else
                        if obs{j}.type==1 % if is pedestrian, then only consider obs avoidance after pedestrian passes
                            if obs{j}.timeL*R*(obs{j}.x0+obs{j}.v*dt*(i-2))+obs{j}.timeL*p+obs{j}.timeS<0 % can start plan to gothrough even have not pass now
%                             if obs{j}.timeL*R*obs{j}.x0+obs{j}.timeL*p+obs{j}.timeS<0 % must wait until pass, then begin to plan go through
                                Lbound=-obs{j}.boundL*R;
                                Sbound=obs{j}.boundL*p+obs{j}.boundS;
                                Lstack=[Lstack;zeros(1,(i-1)*dim) Lbound zeros(1,(nstep-i)*dim) zeros(1,i-1) -1 zeros(1,nstep-i)];
                                Sstack=[Sstack;Sbound-margin];
                            else
                                Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,i-1) -1 zeros(1,nstep-i)];
                                Sstack = [Sstack;S-margin];
                            end
                        else
                            Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,i-1) -1 zeros(1,nstep-i)];
                            Sstack = [Sstack;S-margin];
                        end
                    end
                end
                
%                 Sstack = [Sstack;S];
            end    
            if SCFS==1
                % slack variables
                x2 = path((i-3)*dim+1:(i-2)*dim);
                x1 = path((i-2)*dim+1:(i-1)*dim);
                x0 = path((i-1)*dim+1:i*dim);    
                ur = refinput(i);
%                 % norm^2
%                 du = norm(x1-x2)^2;
%                 dx2 = -2*ur*(x1-x2)+P*(x0-x1);
%                 dx1 = 2*ur*(x1-x2)-P*(x0-x1)+P'*(x1-x2);
%                 dx0 = -P'*(x1-x2);
%                 phi=-(x1-x2)'*P*(x0-x1)+norm(x1-x2)^2*ur;
%                 s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
%                 % norm*norm
%                 du = norm(x1-x2)*norm(x0-x1);
%                 dx2 = -norm(x0-x1)*ur*(x1-x2)/norm(x1-x2)+P*(x0-x1);
%                 dx1 = ur*(-norm(x1-x2)*(x0-x1)/norm(x0-x1)+norm(x0-x1)*(x1-x2)/norm(x1-x2))-P*(x0-x1)+P'*(x1-x2);
%                 dx0 = ur*norm(x1-x2)*(x0-x1)/norm(x0-x1)-P'*(x1-x2);
%                 phi=-(x1-x2)'*P*(x0-x1)+norm(x1-x2)*norm(x0-x1)*ur;
%                 s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
                % curvature
                du = norm(x1-x2)*norm(x0-x1)*norm(x0-x2);
                dx2 = ur*norm(x0-x1)*(norm(x1-x2)*(x2-x0)/norm(x0-x2)+norm(x0-x2)*(x2-x1)/norm(x1-x2))+P*(x0-x1);
                dx1 = ur*norm(x0-x2)*(norm(x1-x2)*(x1-x0)/norm(x0-x1)+norm(x0-x1)*(x1-x2)/norm(x1-x2))-P*(x0-x1)+P'*(x1-x2);
                dx0 = ur*norm(x1-x2)*(norm(x0-x1)*(x0-x2)/norm(x0-x2)+norm(x0-x2)*(x0-x1)/norm(x0-x1))-P'*(x1-x2);
                phi=-(x1-x2)'*P*(x0-x1)+norm(x1-x2)*norm(x0-x1)*norm(x0-x2)*ur;
                s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
                
                Lstack = [Lstack;zeros(1,(i-3)*dim) -dx2' -dx1' -dx0' zeros(1,(nstep-i)*dim) zeros(1,i-1) -du' zeros(1,nstep-i) zeros(1,nstep) zeros(1,nstep)];
                Sstack = [Sstack;s];
                
%                 % norm^2
%                 du = norm(x1-x2)^2;
%                 dx2 = -2*ur*(x1-x2)-P*(x0-x1);
%                 dx1 = 2*ur*(x1-x2)+P*(x0-x1)-P'*(x1-x2);
%                 dx0 = P'*(x1-x2);
%                 phi=(x1-x2)'*P*(x0-x1)+norm(x1-x2)^2*ur;
%                 s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
%                 % norm*norm
%                 du = norm(x1-x2)^2;
%                 dx2 = -norm(x0-x1)*ur*(x1-x2)/norm(x1-x2)-P*(x0-x1);
%                 dx1 = ur*(-norm(x1-x2)*(x0-x1)/norm(x0-x1)+norm(x0-x1)*(x1-x2)/norm(x1-x2))+P*(x0-x1)-P'*(x1-x2);
%                 dx0 = ur*norm(x1-x2)*(x0-x1)/norm(x0-x1)+P'*(x1-x2);
%                 phi=(x1-x2)'*P*(x0-x1)+norm(x1-x2)*norm(x0-x1)*ur;
%                 s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
                % curvature
                du = norm(x1-x2)*norm(x0-x1)*norm(x0-x2);
                dx2 = ur*norm(x0-x1)*(norm(x1-x2)*(x2-x0)/norm(x0-x2)+norm(x0-x2)*(x2-x1)/norm(x1-x2))-P*(x0-x1);
                dx1 = ur*norm(x0-x2)*(norm(x1-x2)*(x1-x0)/norm(x0-x1)+norm(x0-x1)*(x1-x2)/norm(x1-x2))+P*(x0-x1)-P'*(x1-x2);
                dx0 = ur*norm(x1-x2)*(norm(x0-x1)*(x0-x2)/norm(x0-x2)+norm(x0-x2)*(x0-x1)/norm(x0-x1))+P'*(x1-x2);
                phi=(x1-x2)'*P*(x0-x1)+norm(x1-x2)*norm(x0-x1)*norm(x0-x2)*ur;
                s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
                
                Lstack = [Lstack;zeros(1,(i-3)*dim) -dx2' -dx1' -dx0' zeros(1,(nstep-i)*dim) zeros(1,i-1) -du' zeros(1,nstep-i) zeros(1,nstep) zeros(1,nstep)];
                Sstack = [Sstack;s];
                
                % for cos constraint
                du = 2*norm(x0-x1)*norm(x1-x2);
                dx2 = 2*ur*norm(x0-x1)*(x2-x1)/norm(x1-x2)...
                    -2*(x2-x0)+2*(x2-x1)+norm(x0-x1)*(x2-x1)/norm(x1-x2);
                dx1 = 2*ur*(norm(x0-x1)*(x1-x2)/norm(x1-x2)+norm(x1-x2)*(x1-x0)/norm(x0-x1))...
                    +2*(x1-x0)+2*(x1-x2)+(norm(x0-x1)*(x1-x2)/norm(x1-x2)+norm(x1-x2)*(x1-x0)/norm(x0-x1));
                dx0 = 2*ur*norm(x1-x2)*(x0-x1)/norm(x0-x1)...
                    -2*(x0-x2)+2*(x0-x1)+norm(x1-x2)*(x0-x1)/norm(x0-x1);
                phi=-norm(x0-x2)^2+norm(x0-x1)^2+norm(x1-x2)^2+norm(x0-x1)*norm(x1-x2)+2*norm(x0-x1)*norm(x1-x2)*ur;
                s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
                
                Lstack = [Lstack;zeros(1,(i-3)*dim) -dx2' -dx1' -dx0' zeros(1,(nstep-i)*dim) zeros(1,nstep) zeros(1,i-1) -du' zeros(1,nstep-i) zeros(1,nstep)];
                Sstack = [Sstack;s];
                
                du = 2*norm(x0-x1)*norm(x1-x2);
                dx2 = 2*ur*norm(x0-x1)*(x2-x1)/norm(x1-x2)...
                    +2*(x2-x0)-2*(x2-x1)-norm(x0-x1)*(x2-x1)/norm(x1-x2);
                dx1 = 2*ur*(norm(x0-x1)*(x1-x2)/norm(x1-x2)+norm(x1-x2)*(x1-x0)/norm(x0-x1))...
                    -2*(x1-x0)-2*(x1-x2)-(norm(x0-x1)*(x1-x2)/norm(x1-x2)+norm(x1-x2)*(x1-x0)/norm(x0-x1));
                dx0 = 2*ur*norm(x1-x2)*(x0-x1)/norm(x0-x1)...
                    +2*(x0-x2)-2*(x0-x1)-norm(x1-x2)*(x0-x1)/norm(x0-x1);
                phi=norm(x0-x2)^2-norm(x0-x1)^2-norm(x1-x2)^2-norm(x0-x1)*norm(x1-x2)+2*norm(x0-x1)*norm(x1-x2)*ur;
                s = phi-dx2'*x2-dx1'*x1-dx0'*x0-du'*ur;
                
                Lstack = [Lstack;zeros(1,(i-3)*dim) -dx2' -dx1' -dx0' zeros(1,(nstep-i)*dim) zeros(1,nstep) zeros(1,i-1) -du' zeros(1,nstep-i) zeros(1,nstep)];
                Sstack = [Sstack;s];
            end
        end
    end
    options=optimset('Display','off');
    % QP
    if SCFS==1
        % u_sin
        Lstack=[Lstack;zeros(nstep,2*nstep) eye(nstep) zeros(nstep) zeros(nstep)];
%         Sstack=[Sstack;ubar*ones(nstep,1)]; % sin(uq)<=ubar
        Sstack=[Sstack;ubar/2*ones(nstep,1)]; % sin(uq)/v<=ubar*dt
        % u_cos
        Lstack=[Lstack;zeros(nstep,2*nstep) zeros(nstep) 0*eye(nstep) zeros(nstep)];
        Sstack=[Sstack;0.5*ones(nstep,1)];
        % s
        Lstack=[Lstack;zeros(nstep,4*nstep) -eye(nstep)]; %x,u_sin,u_cos,s
        Sstack=[Sstack;zeros(nstep,1)];
%         % bound
%         S=-boundL*p-boundS;L=boundL*R;
%         Lstack=[Lstack;kron(eye(nstep),L) zeros(nstep,3*nstep)];%x,u_sin,u_cos,s
%         Sstack=[Sstack;S*ones(nstep,1)];
        Qe = blkdiag(Qref+Qabs,c(3)*eye(nstep),zeros(nstep),c(4)*eye(nstep));
        f=[-Qref*refpath(:);zeros(2*nstep,1)];
        soln = quadprog(Qe,[f;zeros(nstep,1)],Lstack,Sstack,Aeq,beq,[],[],[],options);
        pathnew = [soln(1:dim*nstep);soln(end-nstep+1:end)];
        refinput = soln(dim*nstep+1:end-nstep);
    else
        Lstack=[Lstack;zeros(nstep,2*nstep) -eye(nstep)];
        Sstack=[Sstack;zeros(nstep,1)];
%         % bound
%         S=-boundL*p-boundS;L=boundL*R;
%         Lstack=[Lstack;kron(eye(nstep),L) zeros(nstep)];
%         Sstack=[Sstack;S*ones(nstep,1)];
        Qe = blkdiag(Qref+Qabs,c(4)*eye(nstep));
        f=-Qref*refpath(:);
        soln = quadprog(Qe,[f;zeros(nstep,1)],Lstack,Sstack,Aeq,beq,[],[],[],options);
%         pathnew = [soln(1:dim*nstep);soln(end-nstep+1:end)];
        pathnew = soln;
    end
    if norm(path-pathnew)/ norm(path)< convergeThres % only based on the 2D path
%         disp(['converged at step ',num2str(k)]);
        break     
    end
    path = pathnew;
end
traj=reshape(path(1:2*nstep),2,nstep);
ss=path(end-nstep+1:end);

end

