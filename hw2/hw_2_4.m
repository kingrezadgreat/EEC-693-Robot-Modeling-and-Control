clc;clear all; close all;
format short 
format compact

%% computing rotations
syms q1 q2 q3 q4 q5
H10 = Rotz(pi/2)*Transz(q1)*Transx(0)*Rotx(pi/2);
H21 = Rotz(pi/2)*Transz(q2)*Transx(0)*Rotx(pi/2);
H32 = Rotz(q3)*Transz(0)*Transx(0)*Rotx(pi/2);
H43 = Rotz(q4)*Transz(0)*Transx(0)*Rotx(-pi/2);
H54 = Rotz(q5)*Transz(0)*Transx(1)*Rotx(0);

H50 = H10*H21*H32*H43*H54;
x(1) = 0;
y(1) = 0;
z(1) = 0;

i = 1;
for t=0:.01:20
    q1 = sin(2*t);
    q2 = cos(t);
    q3 = t;
    q4 = sin(2*t);
    q5 = t;
    H = eval(H50);
    p55 = [0,0,0,1];
    p05 = (H*p55')';
    x(i) = p05(1,1);
    y(i) = p05(1,2);
    z(i) = p05(1,3);
    time(i)=t;
    i = i+1;
    t
end


%% PART 2.5
subplot(2,2,1);
plot(x,time,'.')
xlabel('Time(s)')
ylabel('X')
title('X-time plot')
grid on

subplot(2,2,2);
plot(y,time,'.')
xlabel('Time(s)')
ylabel('Y')
title('Y-time plot')
grid on

subplot(2,2,3);
plot(z,time,'.')
xlabel('Time(s)')
ylabel('Z')
title('Z-time plot')
grid on

subplot(2,2,4);
plot3(x,y,z,'.')
xlabel('X')
ylabel('Y')
zlabel('Z')
title('XYZ plot')
grid on


%% PART 2.6
figure
plot3(x,y,z,'.')
xlabel('X')
ylabel('Y')
zlabel('Z')
title('XYZ plot')
grid on


%% PART 2.7
figure
subplot(2,2,1);
plot(x,y,'.')
xlabel('X')
ylabel('Y')
title('XY projection')
grid on

subplot(2,2,2);
plot(x,z,'.')
xlabel('X')
ylabel('Z')
title('XZ projection')
grid on

subplot(2,2,3);
plot(y,z,'.')
xlabel('Y')
ylabel('Z')
title('YZ projection')
grid on

subplot(2,2,4);
plot3(x,y,z,'.')
xlabel('X')
ylabel('Y')
zlabel('Z')
title('XYZ plot')
grid on
