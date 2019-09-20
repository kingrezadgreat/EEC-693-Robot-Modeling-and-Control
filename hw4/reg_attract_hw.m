
clc
clear all
close all
% format compact


%% syms
syms x1 x2 p11 p12 p21 p22 g ieq xeq1 k1 k2 k m
V = 0.5*[x1 x2]*[p11 p12; p21 p22]*[x1;x2];
x1dot = x2;
u = ieq - k1*(x1-xeq1)-k2*x2;
x2dot = g-(k/m)*u^2/(2*x1^2);
Vdot = diff(V,x1)*x1dot + diff(V,x2)*x2dot;
simplify (Vdot)


%% evaluation


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

figure
dth = 0.1;
dr = 0.01;
rmax = 0.1;


%% Q1
Q=1000*eye(2);
Q1_eig = eig(Q)
R=1;

K=lqr(A,B,Q,R);

Ac = A-B*K;
Ac1_eig = eig(Ac)
P = lyap(Ac',Q);

x1 = 1;
x2 = 1;
p11 = P(1,1);
p12 = P(1,2);
p21 = P(2,1);
p22 = P(2,2);

k1 = K(1,1);
k2 = K(1,2);


for th=0:dth:2*pi
    flag = 0;
    r = 0;
    while (flag == 0)
        x1 = xeq(1,1) + r* cos(th);
        x2 = xeq(2,1) + r* sin(th);
        Vdot = x2*(p11*x1 + (p12*x2)/2 + (p21*x2)/2) + (g - (k*(k2*x2 - ieq + k1*(x1 - xeq1))^2)/(2*m*x1^2))*((p12*x1)/2 + (p21*x1)/2 + p22*x2);
        if Vdot>00.001 || r>rmax
            flag = 1;
        elseif Vdot<0
            plot(x1,x2,'xk')
            hold on
            r = r + dr;
        end
    end
end


%% Q2
% Q=500*eye(2);
Q = [10,2;2,10];
Q2_eig = eig(Q)
R=1;

K=lqr(A,B,Q,R);

Ac = A-B*K;
Ac1_eig = eig(Ac)
P = lyap(Ac',Q);

x1 = 1;
x2 = 1;
p11 = P(1,1);
p12 = P(1,2);
p21 = P(2,1);
p22 = P(2,2);

k1 = K(1,1);
k2 = K(1,2);


for th=0:dth:2*pi
    flag = 0;
    r = 0;
    while (flag == 0 && r<rmax)
        x1 = xeq(1,1) + r* cos(th);
        x2 = xeq(2,1) + r* sin(th);
        Vdot = x2*(p11*x1 + (p12*x2)/2 + (p21*x2)/2) + (g - (k*(k2*x2 - ieq + k1*(x1 - xeq1))^2)/(2*m*x1^2))*((p12*x1)/2 + (p21*x1)/2 + p22*x2);
        if Vdot>0
            flag = 1;
        elseif Vdot<0
            plot(x1,x2,'+r')
            hold on
            r = r + dr;
        end
    end
end



%%
xlabel('x1')
ylabel('x2')
title('Black=Q1 / Red=Q2')
grid on



