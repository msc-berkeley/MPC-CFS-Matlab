%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: derivative
% Function: compute the value and gradient of the 
% Lp signed distance
% Creator: Jianyu Chen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ d, grad] = derivative( x,l )
extend=100;
l(:,1)=l(:,1)+(l(:,1)-l(:,2))/norm(l(:,1)-l(:,2))*extend;
l(:,end)=l(:,end)+(l(:,end)-l(:,end-1))/norm(l(:,end)-l(:,end-1))*extend;
N=size(l,2);
c=3/4;
nx=size(x,2);
grad=zeros(2,nx);
d=zeros(1,nx);
for j=1:nx
    flag=0;
    phi=0;
    dphi=[0;0];
    for i=1:N-1
        
        avec=l(:,i)-x(:,j);
        bvec=l(:,i+1)-x(:,j);
        a=norm(avec);
        b=norm(bvec);
        if a==0||b==0 % when it's on the vertex
            if a==0
                v1=(l(:,i)-l(:,i-1))/norm(l(:,i)-l(:,i-1));
                v2=(l(:,i+1)-l(:,i))/norm(l(:,i+1)-l(:,i));
                v1=[-v1(2);v1(1)];
                v2=[-v2(2);v2(1)];
                grad(:,j)=(v1+v2)/norm(v1+v2);
            elseif b==0
                v1=(l(:,i+1)-l(:,i))/norm(l(:,i+1)-l(:,i));
                v2=(l(:,i+2)-l(:,i+1))/norm(l(:,i+2)-l(:,i+1));
                v1=[-v1(2);v1(1)];
                v2=[-v2(2);v2(1)];
                grad(:,j)=(v1+v2)/norm(v1+v2);            
            end
%             return;
            flag=1;
            break;
        else
            if (avec/a)'*(bvec/b)==-1 % when it's on the line
                grad(:,j)=(l(:,i+1)-l(:,i))/norm(l(:,i+1)-l(:,i));
                grad(:,j)=[-grad(2,j);grad(1,j)]; % rotate ab counter-clockwise for 90 degree           
%                 return;
                flag=1;
                break;
            else % when neither on the line nor on the vertex
                t=tan(acos((bvec/b)'*(avec/a))/2);
                da=-avec/a;db=-bvec/b;
                dt=1/(cos(acos(da'*db)/2))^2*(-1/(2*sqrt(1-(da'*db)^2)))*...
                    (-(avec+bvec)/(a*b)+avec'*bvec*(avec/(a^3*b)+bvec/(b^3*a)));
                sym=cross([avec;0],[bvec;0]);
                if sym(3)>0
                    phi=phi+t/3*(1/a^3+1/b^3)+(t+t^3)/6*(1/a+1/b)^3;
                    if da'*db~=1
                        dphi=dphi+(-t*(da/a^4+db/b^4)+dt/3*(1/a^3+1/b^3)+...
                            (dt+3*t^2*dt)/6*(1/a+1/b)^3-(t+t^3)/2*(1/a+1/b)^2*(da/a^2+db/b^2));
                    end
                else
                    phi=phi-(t/3*(1/a^3+1/b^3)+(t+t^3)/6*(1/a+1/b)^3);
                    if da'*db~=1
                        dphi=dphi-(-t*(da/a^4+db/b^4)+dt/3*(1/a^3+1/b^3)+...
                            (dt+3*t^2*dt)/6*(1/a+1/b)^3-(t+t^3)/2*(1/a+1/b)^2*(da/a^2+db/b^2));     
                    end
                end
            end
        end
    end
    if flag==1
        continue;
    end
    if phi>0
        d(j)=(c*abs(phi))^(-1/3);
    else
        d(j)=-(c*abs(phi))^(-1/3);
    end
    
    grad(:,j)=-1/3*c*d(j)^4*dphi;
end
grad=normc(real(grad));

end