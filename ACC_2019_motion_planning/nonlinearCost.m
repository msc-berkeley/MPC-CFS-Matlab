function [ dL2, dL1 ] = nonlinearCost( x,u,var,ref )
% Add the non-quadratic cost
N=var.N;
nx=var.nx;
nu=var.nu;
vr=var.vr;
dL2=zeros(nx+nu,nx+nu,N);
dL1=zeros(nx+nu,N);
Qcur=zeros(nx+nu);
qcur=zeros(nx+nu,1);
Qref=zeros(nx+nu);
qref=zeros(nx+nu,1);
wcur=var.wcur;
wref=var.wref;
wv=var.wv;
alpha=var.alpha;
% [ a,d ] = findpt( x,ref );
[ d, a] = derivative( x(1:2,:),ref );
% for cost of velocity
as=[0 1;-1 0]*a; 
Qv=zeros(nx+nu);
qv=zeros(nx+nu,1);
for i=1:N-1
    % cost of curvature
    sq=sqrt(u(2,i)^2+alpha^2);
    Qcur(3,3)=2*(sq-alpha);Qcur(3,6)=2*x(3,i)*u(2,i)/sq;Qcur(6,3)=2*x(3,i)*u(2,i)/sq;Qcur(6,6)=(x(3,i)^2*alpha^2)/sq^3;
    qcur(3)=2*x(3,i)*(sq-alpha);qcur(6)=x(3,i)^2*u(2,i)/sq;
    Qcur=Qcur*wcur;qcur=qcur*wcur;
    % cost of ref
    Qref(1:2,1:2)=a(:,i)*a(:,i)';
    qref(1:2)=a(:,i)*d(i);
    Qref=wref*Qref;qref=wref*qref;
    if ~var.pure
        % cost of velocity
        A=as(:,i)*as(:,i)';
        Ad=[2*A(2,1) A(2,2)-A(1,1);A(2,2)-A(1,1) -2*A(2,1)];
        Add=[2*Ad(2,1) Ad(2,2)-Ad(1,1);Ad(2,2)-Ad(1,1) -2*Ad(2,1)];
        B=as(:,i)*vr;
        Bd=[B(2);-B(1)];
        Bdd=[Bd(2);-Bd(1)];
        th=[cos(x(4,i));sin(x(4,i))];
        rvv=th'*A*th; rv=x(3,i)*rvv-th'*B;
        rth=0.5*x(3,i)^2*th'*Ad*th-x(3,i)*th'*Bd;
        rthth=0.5*x(3,i)^2*th'*Add*th-x(3,i)*th'*Bdd;
        rthv=x(3,i)*th'*Ad*th-th'*Bd;
        Qv(3,3)=rvv;Qv(4,4)=rthth;Qv(3,4)=rthv;Qv(4,3)=rthv;qv(3)=rv;qv(4)=rth;
        Qv=wv*Qv;qv=wv*qv;
    else
        % cost of pure velocity tracking
        Qv(3,3)=1;qv(3)=x(3,i)-vr;
        Qv=wv*Qv;qv=wv*qv;
    end
    
    dL2(:,:,i)=Qcur+Qref+Qv;
    dL1(:,i)=qcur+qref+qv;
end
    i=N;
    % cost of ref
    Qref(1:2,1:2)=a(:,i)*a(:,i)';
    qref(1:2)=a(:,i)*d(i);
    Qref=wref*Qref;qref=wref*qref;
    if ~var.pure
        % cost of velocity
        A=as(:,i)*as(:,i)';
        Ad=[2*A(2,1) A(2,2)-A(1,1);A(2,2)-A(1,1) -2*A(2,1)];
        Add=[2*Ad(2,1) Ad(2,2)-Ad(1,1);Ad(2,2)-Ad(1,1) -2*Ad(2,1)];
        B=as(:,i)*vr;
        Bd=[B(2);-B(1)];
        Bdd=[Bd(2);-Bd(1)];
        th=[cos(x(4,i));sin(x(4,i))];
        rvv=th'*A*th; rv=x(3,i)*rvv-th'*B;
        rth=0.5*x(3,i)^2*th'*Ad*th-x(3,i)*th'*Bd;
        rthth=0.5*x(3,i)^2*th'*Add*th-x(3,i)*th'*Bdd;
        rthv=x(3,i)*th'*Ad*th-th'*Bd;
        Qv(3,3)=rvv;Qv(4,4)=rthth;Qv(3,4)=rthv;Qv(4,3)=rthv;qv(3)=rv;qv(4)=rth;
        Qv=wv*Qv;qv=wv*qv;
    else
        % cost of pure velocity tracking
        Qv(3,3)=1;qv(3)=x(3,N)-vr;
        Qv=wv*Qv;qv=wv*qv;
    end

    dL2(:,:,i)=Qcur+Qref+Qv;
    dL1(:,i)=qcur+qref+qv;
    

end

