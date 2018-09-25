function [fun,grad] = cost_curv(x,refpath,Qref,Qabs,Qcurv,nstep)
l = zeros(nstep,1);
dim = size(x,1)/nstep;
P = [0 1;-1 0];
grad = x;
for i=3:nstep
    xb0 = x((i-1)*dim+1:(i-0)*dim);
    xb1 = x((i-2)*dim+1:(i-1)*dim);
    xb2 = x((i-3)*dim+1:(i-2)*dim);
    l(i) = (xb1-xb2)'*P*(xb0-xb1)/norm(xb1-xb2)^2;
end
fun = (x-refpath)'*Qref*(x-refpath)+x'*Qabs*x + l'*Qcurv*l-refpath'*Qref*refpath;
grad = 2*Qref*(x-refpath) + 2*Qabs*x;
end