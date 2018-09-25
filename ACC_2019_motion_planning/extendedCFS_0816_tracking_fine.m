clear all 
close all 

fighandle = [];
fighandle(1) = figure(1); hold on;
set(gcf, 'position', [0 0 500 500]);
fighandle(2) = figure(2); hold on;
fighandle(3) = figure(3); hold on;

%% Set up video
v = VideoWriter('peaks2_st.avi');
v.FrameRate = 5
open(v);

%% Set up state
start = -1;
z0 = [start;0];
zT = [1;0];
vision =1;
size = 0.15;
Vx = 0.15;

trajectory = [];
pathall=[];
pathimplemented=[];

outu = [];
outx = [];

%% Set up obstacles
dt = 0.2;nobj = 3;
var.dt = dt;
obs = {};
obs{1}.v = [0;0];
obs{1}.poly = [0.2 0.8 0.9 0.1; -0.1 -0.1 0.15 0.15];
obs{2}.v = [0;0];
obs{2}.poly = [2 3.2 3.3 1.8; 0.1 0.1 -0.2 -0.2];
obs{3}.v = [0;-0.0];
%obs{3}.poly = [2.8 3.2 3.4 2.9; -0.2 -0.2 0.5 0.5];
obs{3}.poly = [5.3 6.7 6.9 4.8; -0.2 -0.2 0.5 0.5];
move = 3;
figure(fighandle(1));
%plot(path(1,:),path(2,:),'b')
for nobs = 1:nobj-1
    ob = Polyhedron('V',obs{nobs}.poly');
    ob.plot('color','g');
end
for j = 1:nobj        
    poly_now(:,:,j)=obs{j}.poly;
end

axis([-1 8 -2 2]);


%% Set up simulation
nstep = 21;
dim = 2;
SCFS = 1;
ss = 120;
nobjsee = 2;
see  = 15;
%% The cost function
    % The weight
    c = [1,10,20];
    % The distance metric between the original path and the new path
    Q1 = eye(nstep*dim);
    %Q1((nstep-1)*dim+1:end,(nstep-1)*dim+1:end) =  eye(dim)*1000;
    %Q1(1:dim,1:dim) =  eye(dim)*1000;
    % The velocity
    Vdiff = eye(nstep*dim)-diag(ones(1,(nstep-1)*dim),dim);
    Vconst = [-eye(2) eye(2) zeros(2,(nstep-2)*2);[[zeros((nstep-1)*2,2) eye((nstep-1)*2) ]-[eye((nstep-1)*2) zeros((nstep-1)*2,2)]]];
    V_ratio = [5 0;0 1];
    Rpenalty = kron(eye(nstep),V_ratio);
    Q2 = Vconst'*Rpenalty'*Rpenalty*Vconst;
    %Q2 = Vdiff(1:(nstep-1)*dim,:)'*Q1(1+dim:end,1+dim:end)*Vdiff(1:(nstep-1)*dim,:);
    Vref = [Vx,0]*dt;    
    Vref_1 = c(2)*kron(ones(1, nstep),Vref)*Rpenalty'*Rpenalty*Vconst;
    % The accelaration
    Vdiff = eye(nstep*dim)-diag(ones(1,(nstep-1)*dim),dim);
    Adiff = Vdiff-diag(ones(1,(nstep-1)*dim),dim)+diag(ones(1,(nstep-2)*dim),dim*2);
    Q3 = Adiff(1:(nstep-2)*dim,:)'*Adiff(1:(nstep-2)*dim,:);    
   

    
%% MPC
for step = 1:ss
    %% For animation
    %pause(0.2)
    %% Generate reference trajectory
    path = [linspace(z0(1),zT(1),nstep);
            linspace(z0(2),zT(2),nstep)];
    
    
    %% Cost function update
    dir = [zT(1)-(-6) zT(2)-0];
    dir_T = (1/norm(dir))*[ zT(2)-0 -zT(1)+(-6)];
    dd = kron(ones(nstep,1),dir_T*[-6;0]);
    D = kron(eye(nstep),dir_T);
    % Distance to reference line
    Q1 = D'*D;
    Xdis_1 = 2*c(1)*dd'*D;
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
    Aeq = zeros(1*dim,nstep*dim+nstep+nobj*nstep);
    Aeq(0*dim+1:1*dim,1:dim) = eye(dim);      
    beq = [path(:,1)];
    %% The Iteration
    refpath = [];
    for i=1:nstep
        refpath = [refpath;path(:,i)];
    end
    oripath = refpath;
    refinput = ones(1,nstep);
    P = [0 1;-1 0];
    tic
    if step >see
        nobjsee =3
    end
    
%     if step == see-1
%         
%         saveas(figure(1),'2_1_1.eps')
%         pause;
%     elseif step == see+1
%         
%         pause;
%         saveas(figure(1),'2_1_2.eps')
%         pause;
%     elseif step == see+4
%         
%         pause;
%         saveas(figure(1),'2_1_3.eps')
%         pause;
%     end
        
       
    for k = 1:25
    FEAS = 1;
    %% The constraint
    Lstack = []; Sstack = []; margin = size;
   
    for i=1:nstep
        for j=1:nobjsee
            poly = poly_now(:,:,j)+obs{j}.v*ones(1,4)*dt*i;
            [L,S,d] = d2poly(refpath((i-1)*dim+1:i*dim)',poly');
            Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,nstep) zeros(1,(i-1)*nobj+j-1) -1 zeros(1,nobj*(nstep+1-i)-j)];
            Sstack = [Sstack;S-margin]; 
            % Soft constraint
            Lstack = [Lstack;zeros(1,nstep*3) zeros(1,(i-1)*nobj+j-1) -1 zeros(1,nobj*(nstep+1-i)-j)];
            Sstack = [Sstack;0];
        end
        
%         if FEAS > 0 && SCFS > 0
%         if i>2
%             xk0 = refpath((i-1)*dim+1:(i-0)*dim);
%             xk1 = refpath((i-2)*dim+1:(i-1)*dim);
%             xk2 = refpath((i-3)*dim+1:(i-2)*dim);
%             ur = refinput(i);
%             ltheta = norm(xk0-xk1)^2;
%             lk0 = 2*ur*(xk0-xk1)'+(xk1-xk2)'*P';
%             lk1 = -2*ur*(xk0-xk1)'-(xk1-xk2)'*P'+(xk0-xk1)'*P;
%             lk2 = -(xk0-xk1)'*P;
%             s = norm(xk0-xk1)^2*ur+(xk0-xk1)'*P*(xk1-xk2)-ltheta*ur-lk0*xk0-lk1*xk1-lk2*xk2;
%             Lstack = [Lstack;zeros(1,(i-3)*dim) -lk2 -lk1 -lk0 zeros(1,(nstep-i)*dim) zeros(1,i-1) -ltheta zeros(1,nstep-i)];
%             Sstack = [Sstack;s];
%             ltheta = norm(xk0-xk1)^2;
%             lk0 = 2*ur*(xk0-xk1)'-(xk1-xk2)'*P';
%             lk1 = -2*ur*(xk0-xk1)'+(xk1-xk2)'*P'-(xk0-xk1)'*P;
%             lk2 = (xk0-xk1)'*P;
%             s = norm(xk0-xk1)^2*ur-(xk0-xk1)'*P*(xk1-xk2)-ltheta*ur-lk0*xk0-lk1*xk1-lk2*xk2;
%             Lstack = [Lstack;zeros(1,(i-3)*dim) -lk2 -lk1 -lk0 zeros(1,(nstep-i)*dim) zeros(1,i-1) -ltheta zeros(1,nstep-i)];
%             Sstack = [Sstack;s];
%         end
%         end
    end

    %% QP
    if FEAS > 0 && SCFS > 0
        Qe = blkdiag(Qref+Qabs,1*Qcurv);
    else
        Qe = blkdiag(Qref+Qabs,0*Qcurv);
    end
    
    % enlarge A ,f for soft constraint
    Qe = blkdiag(Qe,1000*diag(ones(1,nobj*nstep)));
    
    soln = quadprog(Qe, [-[Xdis_1';zeros(nstep,1)]-[Vref_1';zeros(nstep,1)]; zeros(nobj*nstep,1)],Lstack,Sstack,Aeq,beq);
    pathnew = soln(1:dim*nstep);
    refinput = soln(dim*nstep+1:end);
    z = [pathnew(1*dim+1);pathnew(2*dim)];
    z0 = [pathnew(1*dim+1);pathnew(2*dim)];
    zT = [pathnew(1*dim+1)+vision;0];
    norm(refpath-pathnew)
   
    


    if norm(refpath-pathnew) < 0.1
        pathall = [pathall pathnew]; 
        figure(fighandle(1));
        p(step)=plot(pathnew(1:dim:end),pathnew(2:dim:end),'-*','color',[1-(step/ss),1-(step/ss),1-(step/ss)]);
        axis([-1 8 -2 2]);
        frame = getframe(fighandle(1));
        writeVideo(v,frame);
        
        hold on
        %figure(fighandle(2));
%         plot(refinput,'color',[1-(step/30),1-(step/30),1-(step/30)]);
       
        hold on
        disp(['converged at step ',num2str(k)]);
        
        trajectory = [trajectory z0];
        break
    end
    refpath = pathnew;
    
    end
    
    
    NILQR=20;
    if step ==1
        u0 = zeros(2,NILQR-1);
        x0 = [start;0;0;0;0];
        v0 = Vx;
        angvol = 0;
    else
        u0= zeros(2,NILQR-1);%[u(:,1) zeros(2,NILQR-2)];
        x0 = x(:,2);
        v0 = x(4,2);
        angvol = x(5,2);
    end
    
    traj=zeros(2,NILQR);
    traj(1,:) = pathnew(1:dim:dim*NILQR)';
    traj(2,:) = pathnew(2:dim:dim*NILQR)';
    ILQRtraj = traj(:,1:NILQR);
    xref=[ILQRtraj; zeros(3,NILQR)];
    
   
    
    [ u, x ] = ILQR_tur( x0 , u0, xref, dt );
    z0 = [x(1,2);x(2,2)];
    zT
    zT = [x(1,2)+vision;0]
    
    outu = [outu u(:,1)];
    outx = [outx x(:,2)];
    
    
    plot(z0(1),z0(2),'-o','color','b');
    pp(step) = plot(x(1,:),x(2,:),'-*r');
    pathimplemented = [pathimplemented;z0];
    %%
    time = toc
    disp(['final cost: ']);
    %cost_curv(pathnew,oripath,Qref,Qabs,Qcurv,nstep)  
    
    %% Plot
    ob = Polyhedron('V',obs{3}.poly');        
    if step==see
    PP2 =ob.plot('color','g');
    axis([-1 8 -2 2]);
    end
%     ob = Polyhedron('V',poly_now(:,:,move)');
%     PP =ob.plot('color','b','alpha',((step/(ss)))/2);
%     poly_now(:,:,move) = poly_now(:,:,move)+obs{move}.v*ones(1,4)*dt*step;
     pause(0.1)
    delete(p(step))
    delete(pp(step))
     
%pause;


end
%%

%axis equal
ob = Polyhedron('V',poly_now(:,:,move)');
%PP =ob.plot('color','b','alpha',((step/(ss)))/2);
pend = plot(pathimplemented(1:dim:end),pathimplemented(2:dim:end),'-*','color','r');
%plot(pathall(1:dim:end,10),pathall(2:dim:end,10),'-*','color','g');
grid off
box on
frame = getframe(fighandle(1));
writeVideo(v,frame);
close(v);
%legend([p(1) p(5) p(10) p(20) pend],'Step1','Step5','Step10','Step20','Path implemented','location','best')

% figure(fighandle(2));
% legend('Iter1','Iter2','Iter3','Iter4','Iter5')
% figure;clf; hold on
% curv = [];
% for i=3:nstep
%     l = (refpath((i-2)*dim+1:(i-1)*dim)-refpath((i-3)*dim+1:(i-2)*dim))'*P;
%     s = l*(2*refpath((i-2)*dim+1:(i-1)*dim)-refpath((i-3)*dim+1:(i-2)*dim));
%     v = norm(refpath((i-2)*dim+1:(i-1)*dim)-refpath((i-3)*dim+1:(i-2)*dim))^2;
%     curv(i) = (l*refpath((i-1)*dim+1:i*dim) - s)/v;
% end
% plot(curv);

%% Tragectory
figure
trajectory = pathall';
path = [];
costall =[];
pp = [];
Vconst =[ ];
Rpenalty = [];
Vconst = [-eye(2) eye(2) zeros(2,(ss-2)*2);[[zeros((ss-1)*2,2) eye((ss-1)*2) ]-[eye((ss-1)*2) zeros((ss-1)*2,2)]]];
Rpenalty = kron(eye(ss),V_ratio);
Vdiff = eye(ss*dim)-diag(ones(1,(ss-1)*dim),dim);
Adiff = Vdiff-diag(ones(1,(ss-1)*dim),dim)+diag(ones(1,(ss-2)*dim),dim*2);
Adiff = blkdiag(eye((ss-dim)*dim),zeros(dim*dim))*Adiff  


for i = 1:2:nstep*2
    cost = 0
    path = [];
    path(1:dim:ss*2,1) = trajectory(:,i);
    path(2:dim:ss*2,1) = trajectory(:,i+1);
    I = kron(ones(ss,1),1+(i-1)/2);
    pp(i)=plot3(trajectory(:,i),trajectory(:,i+1),I,'-*','color',[1-(i/(ss*2)),1-(i/(ss*2)),1-(i/(ss*2))]);
    cost = norm(kron(eye(ss),[0 0 ; 0 1])*path)+ norm(Vconst*path-kron(ones(1, ss),Vref)')+norm(Adiff*path)
    costall = [costall  cost]
    hold on
end

zlabel('time')
xlabel('x')
ylabel('y')

figure
plot(costall)
title(' cost VS path ')
xlabel('path_{i}')
ylabel('cost')

%%
k=25;
mm =nstep;
stackpath =[];
diff =[];
mstable =[zeros(mm,1) ones(mm,ss-1)];
for  i = 1:ss
    stackpath(:,i) = [zeros((i-1)*2,1); pathall(:,i);zeros((ss-i)*2,1)];
    
    %error(i-k+M) = norm([trajectory(i,(k-(i-k+M-1))*2-1)-trajectory(i+1,(k-(i-k+M-1))*2-1) trajectory(i,k*2)-trajectory(i+1,k*2)])
end

for i = 1:ss+nstep-1
    for j = 2:ss
        diff(i,j)=norm([stackpath((i-1)*2+1,j)-stackpath((i-1)*2+1,j-1) stackpath(i*2,j)-stackpath(i*2,j-1)]);
    end
end
for M =1:mm
    
    for i = M :ss
        for j =1:M-2
            if diff(i,i-(j-1))<diff(i,i-j) || abs(diff(i,i-(j-1))-diff(i,i-j))<0.005
                mstable(M,i) = mstable(M,i)+1;
            end
        end
    end

%     for i = ss+1 :ss+nstep-1-(M-2)
%         for j =1:M-2
%             if diff(i,20-(j-1))<diff(i,20-j)
%                 mstable(M,i) = mstable(M,i)+1;
%             end
%         end
%     end
end
mstable = max(mstable,[],1)
figure
bar(mstable())
set(gcf, 'position', [0 0 500 200]);
title('M-stable analysis for each action location')
xlabel('k^{th} action location')
ylabel('M-stable')
legend('location','Northwest')
matlab2tikz('myfile.tex');

legend off

figure(fighandle(3));
plot(outx(4,:),'-*g')
plot(outx(5,:),'-*r')
legend('Vol','ANGvol')

figure(4);
plot(outu(1,:),'-*g')
plot(outu(2,:),'-*r')
legend('acc','ANGacc')
    
%legend([pp(1) pp(5) pp(10)  pp(15) pp(20) ],num2str(costall(1)),num2str(costall(5)),num2str(costall(10)),num2str(costall(15)),num2str(costall(20)))