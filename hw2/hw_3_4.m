% Reza Shisheie
% 2708062

clc; 
clear all; 
close all; 
format compact

%% Robot geometry parameters
% please refer to the notes for description of parameters

a = 325;
b = 225;
c = 150;
d = 225;
e = 760;
F = 420+21;
G = F-d-c;
H = 30;
I = 400;
R = 50;
h2 = 100;
h1 = 75;

%% Cylinder geometry onstraints and path planning
xc = I; % x location of the cylinder
yc = 0; % y location of the cylinder
H_w = 40; % height of window
W_w = 15; % width of window
H_b = 10; % offset from base of cylinder
time = 20; % total cut time in seconds
z_l = h1 + H_b; % lower point of the cylinder 
z_h = z_l + H_w; % higher point of cylinder
alpha_max = W_w/R; % maximum cutting angle in radians
L_total = 2*H_w + 2*W_w; % total cut length

V = L_total/time; % linear speed in mm/sec
omega = V/R; % angular speed in rad/sec

alpha_0 = pi/2; % starting angle
z_0 = z_h;

% The assumption is to plot the whole window with a number of dots with constant speed.
% number of dots are 1000
dot_n = 20;
dot_l = floor(dot_n*H_w/L_total);
dot_o = floor(dot_n*W_w/L_total);

i=0;
N1 = 1;
N2 = N1 + dot_o;
N3 = N2 + dot_l;
N4 = N3 + dot_o;
N5 = N4 + dot_l;

z(1) = z_0;
alpha(1) = alpha_0;
tt(1) = 0;
for inc=N1+1:N2
    z(inc) = z_0; 
    alpha(inc) = alpha(inc-1) + alpha_max/ dot_o;
    tt(inc) = tt(inc-1)+time/N5;
    i=i+1;
end

for inc=N2+1:N3
    z(inc) = z(inc-1)-H_w/dot_l; 
    alpha(inc) = alpha(inc-1);
    tt(inc) = tt(inc-1)+time/N5;
    i=i+1;
end

for inc=N3+1:N4
    z(inc) = z(inc-1); 
    alpha(inc) = alpha(inc-1) - alpha_max/ dot_o;
    tt(inc) = tt(inc-1)+time/N5;
    i=i+1;
end

for inc=N4+1:N5
    z(inc) = z(inc-1)+H_w/dot_l; 
    alpha(inc) = alpha(inc-1);
    tt(inc) = tt(inc-1)+time/N5;
    i=i+1;
end

% size(tt)
% size(alpha)
% plot(tt,alpha)
% plot(tt,z)


%% Check one point
% alpha_c = pi/2;
% x_c = I;
% y_c = 0;
% z_c = 175;
% 
% x20 = x_c - (R+H)*cos(alpha_c);
% y20 = y_c - (R+H)*sin(alpha_c);
% 
% K = sqrt(x20^2+y20^2);
% 
% q2 = pi - acos((a^2+b^2-K^2)/(2*a*b));
% q1 = atan2(y20,x20)-atan2((b*sin(q2)),(a+b*cos(q2)));
% q3 = alpha_c - q1 - q2;
% q4 = c + d - z_c;
% 
% 
% q1*180/pi
% q2*180/pi
% q3*180/pi
% q4
% 
% x_check = a*cos(q1) + b*cos(q1+q2) + (H)*cos(q1+q2+q3)
% y_check = a*sin(q1) + b*sin(q1+q2) + (H)*sin(q1+q2+q3)
% z_check = c + d - q4


%% Analytical solution
% inverse kinematics
for i=1:inc
    x20 = xc - (R+H)*cos(alpha(i));
    y20 = yc - (R+H)*sin(alpha(i));
    
    K = sqrt(x20^2+y20^2);
    
    q2(i) = pi - acos((a^2+b^2-K^2)/(2*a*b));
    q1(i) = atan2(y20,x20)-atan2((b*sin(q2(i))),(a+b*cos(q2(i))));
    q3(i) = alpha(i) - q1(i) - q2(i);
    q4(i) = c + d - z(i);
    q_in(:,i) = [q1(i);q2(i);q3(i);q4(i)];
    
end

% check with D-H convention forward kinematics
for i=1:inc
    % check with geometric approach forward kinematics

    x_GA(i) = a*cos(q1(i)) + b*cos(q1(i)+q2(i)) + (H)*cos(q1(i)+q2(i)+q3(i)+pi/2);
    y_GA(i) = a*sin(q1(i)) + b*sin(q1(i)+q2(i)) + (H)*sin(q1(i)+q2(i)+q3(i)+pi/2);
    z_GA(i) = c + d - q4(i);

    % check with D-H convention forward kinematics

    H10 = Rotz(q1(i))*Transz(F)*Transx(a)*Rotx(0);
    H21 = Rotz(q2(i))*Transz(0)*Transx(b)*Rotx(0);
    H32 = Rotz(q3(i))*Transz(-G)*Transx(0)*Rotx(pi);
    H43 = Rotz(0)*Transz(q4(i))*Transx(0)*Rotx(+pi/2);
    H40 = H10*H21*H32*H43;
    p44 = [0,0,H,1];
    p04 = (H40*p44')';
    x_DH(i) = p04(1,1);
    y_DH(i) = p04(1,2);
    z_DH(i) = p04(1,3);
end

figure 
subplot(2,2,1)
plot(x_DH,y_DH,'r*', x_GA,y_GA,'k-.')
xlabel('x(mm)')
ylabel('y(mm)')
grid on
legend('D-H','Geometric')

subplot(2,2,2)
plot(x_DH,z_DH,'r*', x_GA,z_GA,'k-.')
xlabel('x(mm)')
ylabel('z(mm)')
grid on
legend('D-H','Geometric')

subplot(2,2,3)
plot(y_DH,z_DH,'r*', y_GA,z_GA,'k-.')
xlabel('y(mm)')
ylabel('z(mm)')
grid on
legend('D-H','Geometric')

subplot(2,2,4)
plot3(x_DH,y_DH,z_DH,'r*', x_GA,y_GA,z_GA,'k-.')
xlabel('x(mm)')
ylabel('y(mm)')
zlabel('z(mm)')
grid on
%legend('D-H','Geometric')



figure
subplot(2,2,1)
plot(tt,x_DH,'r*', tt,x_GA,'k-.')
xlabel('time(sec)')
ylabel('x_0(mm)')
grid on
legend('D-H','Geometric')

subplot(2,2,2)
plot(tt,y_DH,'r*', tt,y_GA,'k-.')
xlabel('time(sec)')
ylabel('y_0(mm)')
grid on
legend('D-H','Geometric')

subplot(2,2,3)
plot(tt,z_DH,'r*', tt,z_GA,'k-.')
xlabel('time(sec)')
ylabel('z_0(mm)')
grid on
legend('D-H','Geometric')



figure
subplot(2,2,1)
plot(tt,q1,'r')
xlabel('time(sec)')
ylabel('q_1(rad)')
grid on

subplot(2,2,2)
plot(tt,q2,'r')
xlabel('time(sec)')
ylabel('q_2(rad)')
grid on

subplot(2,2,3)
plot(tt,q3,'r')
xlabel('time(sec)')
ylabel('q_3(rad)')
grid on

subplot(2,2,4)
plot(tt,q4,'r')
xlabel('time(sec)')
ylabel('q_4(rad)')
grid on




%% Computation solution

%d1=1;  
%d6=0.5; %location of point of interest on z6 axis

L(1)=Link([0 F a 0]);
L(2)=Link([0 0 b 0]);
L(3)=Link([0 -G 0 pi]);
L(4)=Link([0 0 0 pi/2 1]);
L(5)=Link([0 H 0 0]);

robot=SerialLink(L,'name','myrobot');

W=[0 600 -300 300 0 800]; 
%robot.plot([0 0 0 100 0],'workspace',W); %plot some pose
figure
robot.plot([-pi/4 pi/4 0 100 0],'workspace',W); %plot some pose


%Define Cartesian points for circle using polar parameterization

for i=1:inc
    px(i)=I;
    py(i)=0;
    pz(i)=z(i);

    ox(i)=px(i)-(R)*cos(alpha(i));
    oy(i)=py(i)-(R)*sin(alpha(i));
    oz(i)=pz(i);

end

x0=[1;0;0];y0=[0;1;0];z0=[0;0;1]; %world frame

for i=1:inc
    qguess=[q1(i) q2(i) q3(i) q4(i) 0]; 
    i
    z5=[px(i)-ox(i); py(i)-oy(i);pz(i)-oz(i)]; z5=z5/norm(z5)
    y5=[0;0;-1]; y5=y5/norm(y5)
    x5=cross(y5,z5)
    
    R=[x5'*x0 y5'*x0 z5'*x0;x5'*y0 y5'*y0 z5'*y0;x5'*z0 y5'*z0 z5'*z0];
    H=[R [ox(i);oy(i);oz(i)];0 0 0 1];
    
    qnum(:,i)=robot.ikine(H,qguess,[1 1 1 0 1 0 ])'; 
    qguess=qnum(:,i)'
    
end

% for i=1:inc
%     H=robot.fkine(qnum(:,i)'); %pull out transformation
%     check=H(1:3,4);%look in 4th column to extract world coordinates of endpoint
%     plot3(check(1),check(2),check(3),'ko')
% end
% title('Trajectory: Desired, Analytical and Toolbox')
% xlabel('x_0')
% ylabel('y_0')
% zlabel('z_0')
% 
% 

figure
subplot(2,2,1)
plot(tt,qnum(1,:),'r', tt,q_in(1,:),'k-.')
xlabel('time(sec)')
ylabel('q_1(rad)')
legend('Computation','Analytic')
grid on

subplot(2,2,2)
plot(tt,qnum(2,:),'r', tt,q_in(2,:),'k-.')
xlabel('time(sec)')
ylabel('q_2(rad)')
legend('Computation','Analytic')
grid on

subplot(2,2,3)
plot(tt,qnum(3,:),'r', tt,q_in(3,:),'k-.')
xlabel('time(sec)')
ylabel('q_3(rad)')
legend('Computation','Analytic')
grid on

subplot(2,2,4)
plot(tt,qnum(4,:),'r', tt,q_in(4,:),'k-.')
xlabel('time(sec)')
ylabel('q_4(rad)')
legend('Computation','Analytic')
grid on







