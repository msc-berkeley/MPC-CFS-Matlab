function [] = ellipse( mu, sig )
% plot ellipse from sigma and mu
t=linspace(0,2*pi,1000);
[V,D]=eig(sig);
lamda1=sqrt(D(1,1));lamda2=sqrt(D(2,2));
theta=atan(V(2)/V(1));
x = mu(1) + lamda1*cos(t)*cos(theta) - lamda2*sin(t)*sin(theta);
y = mu(2) + lamda2*sin(t)*cos(theta) + lamda1*cos(t)*sin(theta);
plot(x,y,'r');

end

