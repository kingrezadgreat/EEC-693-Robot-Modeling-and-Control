
clc
clear all
close all

%% step

k = 6.5308*10^(-5);
g = 9.81;
m = 0.068;

xeq1 = 0.006;
xeq2 = 0;
ieq = sqrt(2*m*g/k)*xeq1;

xeq=[xeq1;xeq2];
ueq=ieq;

a = -2*g/ieq;
b = 2*g/xeq1;

A=[0 1;b 0];
B=[0;a];

Q=1000*eye(2);
R=1;

K=lqr(A,B,Q,R);

%Find closed-loop poles
eig(A-B*K)

%Select IC
dx0=[0.001;0];
x0=xeq+dx0;

sim('hw4_step')

figure
subplot(2,2,1)
plot(t,x(:,1)*1000,'r',t,dx(:,1)*1000,'k');
xlabel('time')
ylabel('x1(mm)')
legend('x','dx')
grid on

subplot(2,2,2)
plot(t,x(:,2)*1000,'r',t,dx(:,2)*1000,'k');
xlabel('time')
ylabel('x2(mm/s)')
legend('x','dx')
grid on

subplot(2,2,3)
plot(t,u(:,1),'r',t,du(:,1),'k')
xlabel('time')
ylabel('u')
legend('u','du')
grid on



%% ramp
k = 6.5308*10^(-5);
g = 9.81;
m = 0.068;
xeq1 = 0.014;
xeq2 = 0;
ieq = sqrt(2*m*g/k)*xeq1;

kff = sqrt(2*m*g/k);

xeq=[xeq1;xeq2];
ueq=ieq;

a = -2*g/ieq;
b = 2*g/xeq1;

A=[0 1;b 0];
B=[0;a];

Q=1000*eye(2);
R=1;

K=lqr(A,B,Q,R);

%Find closed-loop poles
eig(A-B*K)

%Select IC
dx0=[0;0];
x0=xeq+dx0;

sim('hw4_ramp')

figure
subplot(3,1,1)
%plot(t,x(:,1)*1000,'r',t,dx(:,1)*1000,'k');
plot(t,x(:,1)*1000);
xlabel('time')
ylabel('x1(mm)')
grid on

subplot(3,1,2)
plot(t,x(:,2)*1000);
xlabel('time')
ylabel('x2(mm/s)')
grid on

subplot(3,1,3)
plot(t,u(:,1))
xlabel('time')
ylabel('u')
grid on
