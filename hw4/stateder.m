function xdot = stateder(x,u)
x1=x(1);
x2=x(2);

k = 6.5308*10^(-5);
g = 9.81;
m = 0.068;

xdot1 = x2;
xdot2 = g - (k/m)*(u^2/(2*x1^2));

% xdot1=-sin(x1)*cos(x2);
% xdot2=-cos(x1)*sin(x2)+3*u;

xdot=[xdot1;xdot2];
