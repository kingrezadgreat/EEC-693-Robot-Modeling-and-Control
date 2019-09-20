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

vpa(H10,2)
vpa(H21,2)
vpa(H32,2)
vpa(H43,2)
vpa(H54,2)

H50 = H10*H21*H32*H43*H54;
vpa(H50,2);

%% PART A
q1 = 1;
q2 = 1;
q3 = pi;
q4 = 0;
q5 = pi/2;

%% PART B
% q1 = 1;
% q2 = -1;
% q3 = pi/2;
% q4 = pi/2;
% q5 = 0;

%% RESULTS
H50 = eval(H50)
p55 = [0,0,0,1]
p05 = (H50*p55')'
