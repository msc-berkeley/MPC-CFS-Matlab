fighandle = [];
fighandle(1) = figure(1); hold on;
fighandle(2) = figure(2); hold on;

path = [-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%  path = [linspace(2.3722,14,N);
%             linspace(-0.9572,0,N)];
%  path = [linspace(3.0928,14,N);
%             linspace(-1.292,0,N)];
    
dt = 1;nobj = 1;
obs = {};
obs{1}.v = [0;0];
obs{1}.poly = [3 6.8 5.8 3.5;3 3 -1 -1];




nstep = size(path,2);
dim = 2;

SCFS = 1;
%% The cost function
% The distance metric between the original path and the new path
Q1 = eye(nstep*dim);
%Q1((nstep-1)*dim+1:end,(nstep-1)*dim+1:end) =  eye(dim)*1000;
%Q1(1:dim,1:dim) =  eye(dim)*1000;
% The velocity
Vdiff = eye(nstep*dim)-diag(ones(1,(nstep-1)*dim),dim);
Q2 = Vdiff(1:(nstep-1)*dim,:)'*Q1(1+dim:end,1+dim:end)*Vdiff(1:(nstep-1)*dim,:);
% The accelaration
Adiff = Vdiff-diag(ones(1,(nstep-1)*dim),dim)+diag(ones(1,(nstep-2)*dim),dim*2);
Q3 = Adiff(1:(nstep-2)*dim,:)'*Adiff(1:(nstep-2)*dim,:);
% The weight
c = [0,2,2];
% The total costj
Qref = 1*(Q1*c(1)+Q2*c(2)+Q3*c(3));
Qabs = 0*Q3*c(3);
%% Extended cost
Mcurv = eye(nstep);
Mcurv(nstep,nstep) = 5;
Vcurv = eye(nstep)-diag(ones(1,nstep-1),1);
Acurv = Vcurv-diag(ones(1,(nstep-1)),1)+diag(ones(1,(nstep-2)),2);
Qcurv = 5*Mcurv;%+Vcurv(1:nstep-1,:)'*Vcurv(1:nstep-1,:)+Acurv(1:(nstep-2),:)'*Acurv(1:(nstep-2),:);
%% The boundary constraint
Aeq = zeros(4*dim,nstep*dim+nstep);
Aeq(0*dim+1:1*dim,1:dim) = eye(dim);
Aeq(1*dim+1:2*dim,(nstep-1)*dim+1:nstep*dim) = eye(dim);
Aeq(2*dim+1:3*dim,1:2*dim) = [-eye(dim) eye(dim)];
Aeq(3*dim+1:4*dim,(nstep-2)*dim+1:nstep*dim) = [-eye(dim) eye(dim)];
beq = [path(:,1);path(:,end);path(:,2)-path(:,1);path(:,end)-path(:,end-1)];
%% The Iteration
refpath = [];
for i=1:nstep
    refpath = [refpath;path(:,i)];
end
oripath = refpath;
refinput = ones(1,nstep);
P = [0 1;-1 0];
tic
for k = 1:10
FEAS = 1;
%% The constraint
Lstack = []; Sstack = []; margin = 0.5;
for i=1:nstep
    for j=1:nobj
        poly = obs{j}.poly+obs{j}.v*ones(1,4)*dt*i;
        [L,S,d] = d2poly(refpath((i-1)*dim+1:i*dim)',poly');
        Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,nstep)];
        Sstack = [Sstack;S-margin];
    end
    if FEAS > 0 && SCFS > 0
    if i>2
        xk0 = refpath((i-1)*dim+1:(i-0)*dim);
        xk1 = refpath((i-2)*dim+1:(i-1)*dim);
        xk2 = refpath((i-3)*dim+1:(i-2)*dim);
        ur = refinput(i);
        ltheta = norm(xk0-xk1)^2;
        lk0 = 2*ur*(xk0-xk1)'+(xk1-xk2)'*P';
        lk1 = -2*ur*(xk0-xk1)'-(xk1-xk2)'*P'+(xk0-xk1)'*P;
        lk2 = -(xk0-xk1)'*P;
        s = norm(xk0-xk1)^2*ur+(xk0-xk1)'*P*(xk1-xk2)-ltheta*ur-lk0*xk0-lk1*xk1-lk2*xk2;
        Lstack = [Lstack;zeros(1,(i-3)*dim) -lk2 -lk1 -lk0 zeros(1,(nstep-i)*dim) zeros(1,i-1) -ltheta zeros(1,nstep-i)];
        Sstack = [Sstack;s];
        ltheta = norm(xk0-xk1)^2;
        lk0 = 2*ur*(xk0-xk1)'-(xk1-xk2)'*P';
        lk1 = -2*ur*(xk0-xk1)'+(xk1-xk2)'*P'-(xk0-xk1)'*P;
        lk2 = (xk0-xk1)'*P;
        s = norm(xk0-xk1)^2*ur-(xk0-xk1)'*P*(xk1-xk2)-ltheta*ur-lk0*xk0-lk1*xk1-lk2*xk2;
        Lstack = [Lstack;zeros(1,(i-3)*dim) -lk2 -lk1 -lk0 zeros(1,(nstep-i)*dim) zeros(1,i-1) -ltheta zeros(1,nstep-i)];
        Sstack = [Sstack;s];
    end
    end
end

%% QP
if FEAS > 0 && SCFS > 0
    Qe = blkdiag(Qref+Qabs,1*Qcurv);
else
    Qe = blkdiag(Qref+Qabs,0*Qcurv);
end
soln = quadprog(Qe,[-Qref*oripath;zeros(nstep,1)],Lstack,Sstack,Aeq,beq);
pathnew = soln(1:dim*nstep);
refinput = soln(dim*nstep+1:end);


figure(fighandle(1));
plot(pathnew(1:dim:end),pathnew(2:dim:end),'-*','color',[1-k/6,1-k/6,1-k/6])
figure(fighandle(2));
plot(refinput,'color',[1-k/6,1-k/6,1-k/6]);
if norm(refpath-pathnew) < 0.1
    disp(['converged at step ',num2str(k)]);
    break
end
refpath = pathnew;
end
%%
time = toc
disp(['final cost: ']);
cost_curv(pathnew,oripath,Qref,Qabs,Qcurv,nstep)

%%
figure(fighandle(1));
%plot(path(1,:),path(2,:),'b')
plot(pathnew(1:dim:end),pathnew(2:dim:end),'k')
ob = Polyhedron('V',obs{1}.poly');
ob.plot('color','g');
axis equal
grid off
box on
legend('Iter1','Iter2','Iter3','Iter4','Iter5')

figure(fighandle(2));
legend('Iter1','Iter2','Iter3','Iter4','Iter5')
% figure;clf; hold on
% curv = [];
% for i=3:nstep
%     l = (refpath((i-2)*dim+1:(i-1)*dim)-refpath((i-3)*dim+1:(i-2)*dim))'*P;
%     s = l*(2*refpath((i-2)*dim+1:(i-1)*dim)-refpath((i-3)*dim+1:(i-2)*dim));
%     v = norm(refpath((i-2)*dim+1:(i-1)*dim)-refpath((i-3)*dim+1:(i-2)*dim))^2;
%     curv(i) = (l*refpath((i-1)*dim+1:i*dim) - s)/v;
% end
% plot(curv);