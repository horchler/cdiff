function cdiff_test
%CDIFF_TEST  

%   Andrew D. Horchler, adh9 @ case . edu, Created 10-30-13
%   Revision: 1.0, 10-30-13


% Using subs
syms t
Phi = [cos(t) cos(t)*sin(t) sin(t)];
Phi_d2 = diff(Phi,t)
double(subs(Phi_d2,t,0))

% Using symfun
syms t
Phi(t) = [cos(t) cos(t)*sin(t) sin(t)];
Phi_d2s = diff(Phi,t)
double(Phi_d2s(0))


% Using complex step diferentiation
Phi = @(t)[cos(t) cos(t).*sin(t) sin(t)];
h = 2^-28;
cdiff = @(f,x)imag(f(x(:)+1i*h))/h;
Phi_d2 = cdiff(Phi,0)

t = linspace(0,2*pi,1000);
Phi_d2 = cdiff(Phi,t);
figure
plot(t',Phi_d2(:,1),'b',t,Phi_d2(:,2),'g',t,Phi_d2(:,3),'r')


t = linspace(0,2*pi,1000);
Phi = [cos(t); cos(t).*sin(t); sin(t)];
dPhi = gradient(Phi,t(2)-t(1));

dPhis = Phi_d2s(t)
dPhis = reshape(double([dPhis{:}]),1000,3)';

figure
plot(t,dPhis(1,:),'b',t,dPhis(2,:),'g',t,dPhis(3,:),'r')
whos

figure
err1 = abs(dPhis-dPhi);
plot(t,err1(1,:),'b',t,err1(2,:),'g',t,err1(3,:),'r')