%% PROBLEM 2
clc; clear all; close all; format short; format compact;

%% PART 1: 
fprintf("\n---PROBLEM2 - PART 1---\n")
% Obtain all transformations necessary to compute the Jacobian. As a verification 
% reference, set q 1 = 1, q 2 = −1, q 3 = π/2, q 4 = π/2, q 5 = 0.

syms q1 q2 q3 q4 q5

H10 = Rotz(pi/2)*Transz(q1)*Transx(0)*Rotx(pi/2);
H21 = Rotz(-pi/2)*Transz(q2)*Transx(0)*Rotx(-pi/2);
H32 = Rotz(q3)*Transz(1)*Transx(0)*Rotx(-pi/2);
H43 = Rotz(q4)*Transz(0)*Transx(0)*Rotx(pi/2);
H54 = Rotz(q5)*Transz(0)*Transx(-1)*Rotx(0);

H20 = H10*H21;
H30 = H20*H32;
H40 = H30*H43;
H50 = H40*H54;

% Computing H50 with initial conditions
q1 = 1; q2 = -1; q3 = pi/2; q4 = pi/2; q5 = 0;
H_eval_check = eval(H50)

%% PART 2: Finding the velocity Jacobian Jv using symbolic processing.
fprintf("\n---PROBLEM2 - PART 2---\n")

% finding center of each coordinate system
o1 = H10(1:3,4);
o2 = H20(1:3,4);
o3 = H30(1:3,4);
o4 = H40(1:3,4);
o5 = H50(1:3,4);

% finding zeros
z0 = [0,0,1];
z1 = H10(1:3,3);
z2 = H20(1:3,3);
z3 = H30(1:3,3);
z4 = H40(1:3,3);
z5 = H50(1:3,3);

% geeting jacobian
Jv1 = z0';
Jv2 = z1;
Jv3 = cross(z2,(o5-o2));
Jv4 = cross(z3,(o5-o3));
Jv5 = cross(z4,(o5-o4));
Jw1 = [0;0;0];
Jw2 = [0;0;0];
Jw3 = z2;
Jw4 = z3;
Jw5 = z4;

Jv = [Jv1, Jv2, Jv3, Jv4, Jv5]; 
Jw = [Jw1, Jw2, Jw3, Jw4, Jw5];
J = [Jv; Jw]

%% PART 3: 
% As a cautionary note, try the rank command on the symbolic J v , with
% unspecified values for q. Then evaluate J v using the test values of q above,
% entering sym(pi) for π. Comment on the results.
fprintf("\n---PROBLEM2 - PART 3---\n")

% here is the rank of Symbolic Jv. It is 3 which means it is full rank as
% far as linear velocity
rank_jv_sym = rank(Jv)

% Here is rank of Jv if it is initialized. Rank is 2 which menas that it has 
% lost one degree of freedom in Jacobian and linear velocity is limited in 
% one direction.
q1 = 1; q2 = -1; q3 = pi/2; q4 = pi/2; q5 = 0;
rank_jv_eval = rank(eval(Jv))

%% PART 4: Based on the above, can J v be singular at certain q values? For 
%  example which ones?

% Jv can be singular at certain qs. One of them is the one mentioned above

%% PART 5: Select a q that produces singularity and provide a physical/geometric
%  explanation of how this occurs.
fprintf("\n---PROBLEM2 - PART 5---\n")

% To have an idea what the physical position of the end effector is lets
% find p once all angles are zero and q1=q2=1
q1 = 1; 
q2 = 1; 
q3 = 0;
q4 = 0; 
q5 = 0;
p = eval(H50*[0 0 0 1]'); 
p_nonsing = p(1:3,1)'
rank_jv_nonsing_eval = rank(eval(Jv))
% p yields [1 1 2] which is an upright position
% this is also a full rank situation

q1 = 1; 
q2 = 1; 
q3 = 0; 
q4 = pi/2; 
q5 = 0;
p = eval(H50*[0 0 0 1]'); 
p_sing = p(1:3,1)'
rank_jv_sing_eval = rank(eval(Jv))
% in thhis orientation p yields [1 2 1]. Rank of Jacobian in this position
% is 2 which is not full rank. q1 and q2 can generate any planar motion in 
% the xz plane. Singularity happens when the rod p is in such a orinetation 
% that the motion of q3 q4 a5 only generate velocity in direction of q1 and
% q2. that is an example of a situation when singularity happens. In this
% example q4 turns 90 degrees such that the rod p is pointing up toward y0
% and in this orientation it can only provide velocity in the xz plane and
% not in y0 direction. This would be a singularity situation where one
% degree of freedom in velocity is lost.

%% PART 6: Compute Yoshikawa’s manipulability measure μ(q) using symbolic 
%  processing. Provide a mathematical and/or physical explanation for the result.
fprintf("\n---PROBLEM2 - PART 6---\n")

% mu_sym = sqrt(det(J.'*J));

% non singular mu
q1 = 1; 
q2 = 1; 
q3 = 0;
q4 = 0; 
q5 = 0;
J_eval = eval(Jv);
mu_nonsing = sqrt(det(J_eval*J_eval.'))
rank_Jv_nosing = rank(J_eval)

% singular mu
q1 = 1; 
q2 = 1; 
q3 = 0; 
q4 = pi/2; 
q5 = 0;
J_eval = eval(Jv);
mu_sing = sqrt(det(J_eval*J_eval.'))
rank_Jv_sing = rank(J_eval)
% as you see the determinant of Yoshikawa manipulability is 0 for singular
% orientation and non-zero for non-singular orientation. The singular and
% non-singular orientation was were already menationed above


%% PART 7: For the above test values of q, describe the set of all joint velocities that
%  result in an instantaneous zero velocity for P . Give an example of such
%  a joint velocity vector (other than the zero vector) and verify it with a
%  forward calculation.

fprintf("\n---PROBLEM2 - PART 7---\n")

% lets compute for a non singular situation. 
q1 = 1; 
q2 = 1; 
q3 = 0; 
q4 = 0; 
q5 = 0;
q1dot = 1; 
q2dot = 2; 
q3dot = 3; 
q4dot = 4; 
q5dot = 5;
p = eval(H50*[0 0 0 1]'); 
p_nonsing = p(1:3,1)'
Jv_eval = eval(Jv);
% forward kinematics
xdot_nonsing = (Jv_eval*[q1dot, q2dot, q3dot, q4dot, q5dot]')'
% as you see all values of xdot is non zero


% lets compute for a singular situation now 
q1 = 1; 
q2 = 1; 
q3 = 0; 
q4 = -pi/2; 
q5 = 0;
q1dot = 1; 
q2dot = 2; 
q3dot = 3; 
q4dot = 4; 
q5dot = 5;
p = eval(H50*[0 0 0 1]'); 
p_sing = p(1:3,1)'
Jv_eval = eval(Jv);
% forward kinematics
xdot_sing = (Jv_eval*[q1dot, q2dot, q3dot, q4dot, q5dot]')'
% The orientation of rod is pointing down and thus the location of p is at
% [1,0,1]. In this orientation rod is in direction of y and our expectation
% is that the velocity in y direction should be zero. As you see all values 
% of xdot is non zero except the y element as we predicted

%% PART 8-9:  
% • Calculate the angular velocity Jacobian J w and exclude the first two
% columns. This is equivalent to considering the wrist only, and only ori-
% entation kinematics.
% • For the above reduced angular velocity Jacobian, compute Yoshikawa’s
% manipulability measure in symbolic form. It predicts singularity for cer-
% tain values of q 4 . Interpret geometrically.
fprintf("\n---PROBLEM2 - PART 8-9---\n")

Jw_wrist = Jw(:,3:5);
q4 = 0;
det_sing_0 = det(eval(Jw_wrist))

Jw_wrist = Jw(:,3:5);
q4 = pi;
det_sing_pi = det(eval(Jw_wrist))

% determinant is zero which implies singularity when 
% q4 = n*pi , n = 0,1,2,3,...
% this is called gimbal lock and this phenomena happens when two rotational
% axis are aligned to each other. When q4 = 0 or pi q3 and q5 are aligned
% and thus one degree of freedom is lost.





















