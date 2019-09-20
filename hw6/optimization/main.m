clc;
clear all;
close all;
format compact;

L1 = 1;
L2 = 10;

X0 = rand(1,24);
t = 0:0.1:10;

%% SET CONSTRAINTS
q1_min = -1;
q1_max = 1;
q2_min = -1;
q2_max = 1;
q4_min = -0.5;
q4_max = 2;

q1d_min = -1;
q1d_max = 1;
q2d_min = -1;
q2d_max = 1;
q4d_min = -1;
q4d_max = 1;

ang_const = [q1_min, q1_max, q2_min, q2_max, q4_min q4_max, q1d_min, q1d_max, q2d_min , q2d_max, q4d_min, q4d_max]; 
% ang_const = ang_const*2;

%% FMINCON WITH CONSTRAINTS
% options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-11,'StepTolerance',1e-11,'MaxIterations',10000);
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-11,'StepTolerance',1e-11,'MaxIterations',10000);
Xsol=fmincon(@(X) cost_function(X,L1,L2,t), X0, [],[], [], [], [],[], @(X) nonlconstr(X,t,ang_const),options);


a1 = Xsol(1:4);
a2 = Xsol(5:8);
a4 = Xsol(9:12);

b1 = Xsol(13:16);
b2 = Xsol(17:20);
b4 = Xsol(21:24);

w = 2*pi/10;
iter = 1;
for tt = t
    
    q1=0;q2=0;q4=0;
    q1d=0;q2d=0;q4d=0;
    q1dd=0;q2dd=0;q4dd=0;
    for j = 1:4
        q1 = q1 + a1(j)*sin(w*j*tt)/(w*j) - b1(j)*cos(w*j*tt)/(w*j);
        q2 = q2 + a2(j)*sin(w*j*tt)/(w*j) - b2(j)*cos(w*j*tt)/(w*j);
        q4 = q4 + a4(j)*sin(w*j*tt)/(w*j) - b4(j)*cos(w*j*tt)/(w*j);
        
        q1d = q1d + a1(j)*cos(w*j*tt) + b1(j)*sin(w*j*tt);
        q2d = q2d + a2(j)*cos(w*j*tt) + b2(j)*sin(w*j*tt);
        q4d = q4d + a4(j)*cos(w*j*tt) + b4(j)*sin(w*j*tt);
        
        q1dd = q1dd - a1(j)*w*j*sin(w*j*tt) - b1(j)*w*j*cos(w*j*tt);
        q2dd = q2dd - a2(j)*w*j*sin(w*j*tt) - b2(j)*w*j*cos(w*j*tt);
        q4dd = q4dd - a4(j)*w*j*sin(w*j*tt) - b4(j)*w*j*cos(w*j*tt);
    end
    
    q(:,iter) = [q1;q2;q4];
    qd(:,iter) = [q1d;q2d;q4d];
    qdd(:,iter) = [q1dd;q2dd;q4dd];
    iter = iter + 1;
end




figure
subplot(3,2,1)
plot(t,q(1,1:length(t)))
ylabel('q1')
subplot(3,2,3)
plot(t,q(2,1:length(t)))
ylabel('q2')
subplot(3,2,5)
plot(t,q(3,1:length(t)))
ylabel('q4')

subplot(3,2,2)
plot(t,qd(1,1:length(t)))
ylabel('q1d')
subplot(3,2,4)
plot(t,qd(2,1:length(t)))
ylabel('q2d')
subplot(3,2,6)
plot(t,qd(3,1:length(t)))
ylabel('q4d')

suptitle('with constraint')


%% FMINCON WITH NO CONSTRAINTS
Xsol=fmincon(@(X) cost_function(X,L1,L2,t), X0, [],[], [], [], [],[], [],options);


a1 = Xsol(1:4);
a2 = Xsol(5:8);
a4 = Xsol(9:12);

b1 = Xsol(13:16);
b2 = Xsol(17:20);
b4 = Xsol(21:24);

w = 2*pi/10;
iter = 1;
for tt = t
    
    q1=0;q2=0;q4=0;
    q1d=0;q2d=0;q4d=0;
    q1dd=0;q2dd=0;q4dd=0;
    for j = 1:4
        q1 = q1 + a1(j)*sin(w*j*tt)/(w*j) - b1(j)*cos(w*j*tt)/(w*j);
        q2 = q2 + a2(j)*sin(w*j*tt)/(w*j) - b2(j)*cos(w*j*tt)/(w*j);
        q4 = q4 + a4(j)*sin(w*j*tt)/(w*j) - b4(j)*cos(w*j*tt)/(w*j);
        
        q1d = q1d + a1(j)*cos(w*j*tt) + b1(j)*sin(w*j*tt);
        q2d = q2d + a2(j)*cos(w*j*tt) + b2(j)*sin(w*j*tt);
        q4d = q4d + a4(j)*cos(w*j*tt) + b4(j)*sin(w*j*tt);
        
        q1dd = q1dd - a1(j)*w*j*sin(w*j*tt) - b1(j)*w*j*cos(w*j*tt);
        q2dd = q2dd - a2(j)*w*j*sin(w*j*tt) - b2(j)*w*j*cos(w*j*tt);
        q4dd = q4dd - a4(j)*w*j*sin(w*j*tt) - b4(j)*w*j*cos(w*j*tt);
    end
    
    q(:,iter) = [q1;q2;q4];
    qd(:,iter) = [q1d;q2d;q4d];
    qdd(:,iter) = [q1dd;q2dd;q4dd];
    iter = iter + 1;
end



figure
subplot(3,2,1)
plot(t,q(1,1:length(t)))
ylabel('q1')
subplot(3,2,3)
plot(t,q(2,1:length(t)))
ylabel('q2')
subplot(3,2,5)
plot(t,q(3,1:length(t)))
ylabel('q4')

subplot(3,2,2)
plot(t,qd(1,1:length(t)))
ylabel('q1d')
subplot(3,2,4)
plot(t,qd(2,1:length(t)))
ylabel('q2d')
subplot(3,2,6)
plot(t,qd(3,1:length(t)))
ylabel('q4d')

suptitle('no constraint')



