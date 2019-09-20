clc;clear all; close all;
format short 
format compact

%% computing rotations
syms q1 q2 q3 q4
R1 = Rx(q1)
R2 = Ry(q2)
R3 = Ry(q3)
R4 = Rx(q4)

%% computing final rotation
Rot = R4*R2*R1*R3

%% evaluating rotation matrix
q1=pi/2
q2=-pi
q3=pi/2
q4=-pi/2
Rot = eval(Rot)
%Res_deg = radtodeg(tr2eul(Rot))

%% computing euler angles
theta1 = atan2( sqrt(1-Rot(3,3)^2), Rot(3,3));
phi1 = atan2(Rot(2,3),Rot(1,3));
psi1 = atan2(Rot(3,2),-Rot(3,1));
res1 = [phi1, theta1, psi1 ]
%res1_d = radtodeg([phi1, phi1, psi1 ])

theta2 = atan2( -sqrt(1-Rot(3,3)^2), Rot(3,3));
phi2 = atan2(-Rot(2,3),-Rot(1,3));
psi2 = atan2(-Rot(3,2),Rot(3,1));
res2 = [phi2, theta2, psi2 ]
%res2_d = radtodeg([phi2, phi2, psi2 ])

%% checking to see if the same rotation is resulted
Rot_check_1 = Rz(phi1)*Ry(theta1)*Rz(psi1) - Rot
Rot_check_2 = Rz(phi2)*Ry(theta2)*Rz(psi2) - Rot

