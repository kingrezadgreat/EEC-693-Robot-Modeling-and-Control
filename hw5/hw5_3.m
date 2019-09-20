%Kinematics and parameters for the WAM 4 DOF robot
%Third joint will be frozen at q3=0

%Frames and kinematics:
clc,
clear all;
close all;
format compact

%% 1,2. Kinematics
syms q1 q2 q3 q4 a d

H10=Rotz(q1)*Rotx(-sym(pi)/2);
H21=Rotz(q2)*Rotx(sym(pi)/2);
H32=Rotz(q3)*Transz(d)*Transx(a)*Rotx(-sym(pi)/2);
H43=Rotz(q4)*Transx(-a)*Rotx(sym(pi)/2);  

H20=H10*H21;
H30=H20*H32;
H40=H30*H43;

%Rotations
R1=H10(1:3,1:3);
R2=H20(1:3,1:3);
R3=H30(1:3,1:3);
R4=H40(1:3,1:3);

%z axes
z0=[0;0;1];
z1=H10(1:3,3);
z2=H20(1:3,3);
z3=H30(1:3,3);
z4=H40(1:3,3);

%CM locations in local frames
%cmxx denotes cm of link x in frame x]

syms xc1 yc1 zc1 xc2 zc2 yc2 xc3 yc3 zc3 xc4 yc4 zc4

cm11=[xc1;yc1;zc1;1];
cm22=[xc2;yc2;zc2;1];
cm33=[xc3;yc3;zc3;1];
cm44=[xc4;yc4;zc4;1];

%CM locations in world frame
cm10=H10*cm11;cm10=cm10(1:3);
cm20=H20*cm22;cm20=cm20(1:3);
cm30=H30*cm33;cm30=cm30(1:3);
cm40=H40*cm44;cm40=cm40(1:3);

%Origin locations in world frame
o10=H10(1:3,4);
o20=H20(1:3,4);
o30=H30(1:3,4);
o40=H40(1:3,4);

%Velocity Jacobians relative to CMs
%Link 1
J1v1=cross(z0,cm10);
J1v=[J1v1 zeros(3,3)];

%Link 2
J2v1=cross(z0,cm20);
J2v2=cross(z1,cm20-o10);
J2v=[J2v1 J2v2 zeros(3,2)];

%Link 3
J3v1=cross(z0,cm30);
J3v2=cross(z1,cm30-o10);
J3v3=cross(z2,cm30-o20);
J3v=[J3v1 J3v2 J3v3 zeros(3,1)];

%Link 4
J4v1=cross(z0,cm40);
J4v2=cross(z1,cm40-o10);
J4v3=cross(z2,cm40-o20);
J4v4=cross(z3,cm40-o30);
J4v=[J4v1 J4v2 J4v3 J4v4];

%Angular velocity Jacobians
J1w=[z0 zeros(3,3)];
J2w=[z0 z1 zeros(3,2)];
J3w=[z0 z1 z2 zeros(3,1)];
J4w=[z0 z1 z2 z3];

%% 3. Obtain Mass matrix

syms xc1 yc1 zc1 xc2 yc2 zc2 xc3 yc3 zc3 xc4 yc4 zc4
syms I1xx I1xy I1xz I1yy I1yz I1zz
syms I2xx I2xy I2xz I2yy I2yz I2zz
syms I3xx I3xy I3xz I3yy I3yz I3zz
syms I4xx I4xy I4xz I4yy I4yz I4zz

syms m1 m2 m3 m4

I1 = [I1xx, I1xy, I1xz;
      I1xy, I1yy, I1yz;
      I1xz, I1yz, I1zz];

I2 = [I2xx, I2xy, I2xz;
      I2xy, I2yy, I2yz;
      I2xz, I2yz, I2zz];
  
  
I3 = [I3xx, I3xy, I3xz;
      I3xy, I3yy, I3yz;
      I3xz, I3yz, I3zz];

I4 = [I4xx, I4xy, I4xz;
      I4xy, I4yy, I4yz;
      I4xz, I4yz, I4zz];
  
  

sum = J1w.'*R1*I1*R1.'*J1w + m1*J1v.'*J1v;
sum = sum + m2*J2v.'*J2v + J2w.'*R2*I2*R2.'*J2w;
sum = sum + m3*J3v.'*J3v + J3w.'*R3*I3*R3.'*J3w;
sum = sum + m4*J4v.'*J4v + J4w.'*R4*I4*R4.'*J4w;

% sum = simplify(sum);
% sum(1,1)

%% 4. Set q 3 = 0 and evaluate M to take in this value. Then remove the 
% third row and the third column to obtain a 3-by-3 matrix.
q3 = 0;
sum = eval(sum);
M = sum;
D_q=simplify(sum);
D_q(3,:) = [];
D_q(:,3) = [];



%% 5. Obtain the Coriolis matrix symbolically using M and the gravity vector.
q=[q1;q2;q4];
for i=1:3
    for j=1:3
        for k=1:3
            c(i,j,k)=(1/2)*((diff(D_q(k,j),q(i)))+(diff(D_q(k,i),q(j)))-(diff(D_q(i,j),q(k))));
        end
    end
end

syms q1dot q2dot q4dot
qdot = [q1dot;q2dot;q4dot];

for k=1:3
    for j=1:3
        C(k,j)=c(1,j,k)*q1dot+c(2,j,k)*q2dot+c(3,j,k)*q4dot;
    end
end

%The above has been checked to satisfy the skew-symmetry property for
%Ddot-2C

%Potential Energy
syms g
gv=[0;0;g];  %gravity direction in world frame
P=-m1*gv.'*cm10 - m2*gv.'*cm20 - m3*gv.'*cm30 - m4*gv.'*cm40;

gq1=-simplify(diff(P,q1));
gq2=-simplify(diff(P,q2));
gq3=-simplify(diff(P,q3));
gq4=-simplify(diff(P,q4));


gq=[gq1;gq2;gq3;gq4];
q3 = 0;
gq(3,:) = [];
gq=eval(gq);



%% 6. Verify the skew-symmetry property by symbolic differentiation of M
D_q_skew_symmetry = simplify(D_q.'-D_q)


%% 7.Inspect each one of the 6 distinct entries of M to identify Θ parameters.

syms w
yy(1,1) = w;
D_q(1,1)

TH1 = I2xx/2 + I3xx/2 + I4xx/2 + I1yy + I3yy/2 + I2zz/2 + I4zz/2;
% D_q(1,1) = TH1 - D_q(1,1);
yy(1,1) = 1;
term1 = D_q(1,1) - TH1

yy(1,2) = cos(2*q2 + 2*q4);
term =  subs(term1,yy(1,2),w);
TH2 = diff(term,w);
term2 = simplify(term-TH2*w)

yy(1,3) = sin(2*q2 + 2*q4);
term =  subs(term2,yy(1,3),w);
TH3 = diff(term,w);
term3 = simplify(term-TH3*w)

yy(1,4) = cos(2*q2);
term =  subs(term3, yy(1,4),w);
TH4 = diff(term,w);
term4 = simplify(term-TH4*w)

yy(1,5) = sin(2*q2);
term =  subs(term4, yy(1,5), w);
TH5 = diff(term,w);
term5 = simplify(term-TH5*w)

yy(1,6) = cos(q4);
term =  subs(term5, yy(1,6), w);
TH6 = diff(term,w);
term6 = simplify(term-TH6*w)

yy(1,7) = cos(2*q2 + q4);
term =  subs(term6, yy(1,7), w);
TH7 = diff(term,w);
term7 = simplify(term-TH7*w)

yy(1,8) = sin(2*q2 + q4);
term =  subs(term7, yy(1,8), w);
TH8 = diff(term,w);
term8 = simplify(term-TH8*w)

yy(1,9) = sin(q4);
term =  subs(term8, yy(1,9), w);
TH9 = diff(term,w);
term9 = simplify(term-TH9*w)

yy(1,10) = 1;
TH10 = (a^2*m3)/2 + a^2*m4 + (d^2*m3)/2 + (d^2*m4)/2 + m1*xc1^2 + (m2*xc2^2)/2 + (m3*xc3^2)/2 + (m4*xc4^2)/2 + m2*yc2^2 + (m3*yc3^2)/2 + m4*yc4^2 + m1*zc1^2 + (m2*zc2^2)/2 + m3*zc3^2 + (m4*zc4^2)/2 + a*m3*xc3 - a*m4*xc4 - d*m3*yc3;
TH(10,1) = TH10; 
term10 = term9-TH10

D_11_TH = (TH1+TH2*yy(1,2)+TH3*yy(1,3)+TH4*yy(1,4)+TH5*yy(1,5)+TH6*yy(1,6)+TH7*yy(1,7)+TH8*yy(1,8)+TH9*yy(1,9)+TH10);
diff_check_1 = eval(D_q(1,1)-D_11_TH);
diff_check_1 = simplify(diff_check_1)

% D_q(1,2)
term10 = D_q(1,2)
yy(1,11) = cos(q2 + q4);
term =  subs(term10, yy(1,11), w);
TH11 = diff(term,w);
term11 = simplify(term-TH11*w)

yy(1,12) = cos(q2);
term =  subs(term11, yy(1,12), w);
TH12 = diff(term,w);
term12 = simplify(term-TH12*w)

yy(1,13) = sin(q2 + q4);
term =  subs(term12, yy(1,13), w);
TH13 = diff(term,w);
term13 = simplify(term-TH13*w)

yy(1,14) = sin(q2);
term =  subs(term13, yy(1,14), w);
TH14 = diff(term,w);
term14 = simplify(term-TH14*w)

D_12_TH = (TH11*yy(1,11)+TH12*yy(1,12)+TH13*yy(1,13)+TH14*yy(1,14));
diff_check_2 = eval(D_q(1,2)-D_12_TH);
diff_check_2 = simplify(diff_check_2)



% D_q(1,3)
term14 = D_q(1,3)
yy(1,15) = cos(q2 + q4);
term =  subs(term14, yy(1,15), w);
TH15 = diff(term,w);
term15 = simplify(term-TH15*w)

yy(1,16) = sin(q2 + q4);
term =  subs(term15, yy(1,16), w);
TH16 = diff(term,w);
term16 = simplify(term-TH16*w)

D_13_TH = (TH15*yy(1,15)+TH16*yy(1,16));
diff_check_3 = eval(D_q(1,3)-D_13_TH);
diff_check_3 = simplify(diff_check_3)


% D_q(2,1)
term16 = D_q(2,1)
yy(1,17) = cos(q2 + q4);
term =  subs(term16, yy(1,17), w);
TH17 = diff(term,w);
term17 = simplify(term-TH17*w)

yy(1,18) = cos(q2);
term =  subs(term17, yy(1,18), w);
TH18 = diff(term,w);
term18 = simplify(term-TH18*w)

yy(1,19) = sin(q2 + q4);
term =  subs(term18, yy(1,19), w);
TH19 = diff(term,w);
term19 = simplify(term-TH19*w)

yy(1,20) = sin(q2);
term =  subs(term19, yy(1,20), w);
TH20 = diff(term,w);
term20 = simplify(term-TH20*w)

D_21_TH = (TH17*yy(1,17)+TH18*yy(1,18)+TH19*yy(1,19)+TH20*yy(1,20));
diff_check_4 = eval(D_q(2,1)-D_21_TH);
diff_check_4 = simplify(diff_check_4)


% D_q(2,2)
term20 = D_q(2,2)
yy(1,21) = cos(q4);
term =  subs(term20, yy(1,21), w);
TH21 = diff(term,w);
term21 = simplify(term-TH21*w)

yy(1,22) = sin(q4);
term =  subs(term21, yy(1,22), w);
TH22 = diff(term,w);
term22 = simplify(term-TH22*w)

yy(1,23) = 1;
TH23 = I2yy + I4yy + I3zz + a^2*m3 + 2*a^2*m4 + d^2*m3 + d^2*m4 + m2*xc2^2 + m3*xc3^2 + m4*xc4^2 + m3*yc3^2 + m2*zc2^2 + m4*zc4^2 + 2*a*m3*xc3 - 2*a*m4*xc4 - 2*d*m3*yc3;
term23 = term22-TH23

D_22_TH = (TH21*yy(1,21)+TH22*yy(1,22)+TH23);
diff_check_5 = eval(D_q(2,2)-D_22_TH);
diff_check_5 = simplify(diff_check_5)



% D_q(2,3)
term23 = D_q(2,3)
yy(1,24) = cos(q4);
term =  subs(term23, yy(1,24), w);
TH24 = diff(term,w);
term24 = simplify(term-TH24*w)

yy(1,25) = sin(q4);
term =  subs(term24, yy(1,25), w);
TH25 = diff(term,w);
term25 = simplify(term-TH25*w)

yy(1,26) = 1;
TH26 = m4*a^2 - 2*m4*a*xc4 + m4*xc4^2 + m4*zc4^2 + I4yy;
term26 = term25-TH26

D_23_TH = (TH24*yy(1,24)+TH25*yy(1,25)+TH26);
diff_check_6 = eval(D_q(2,3)-D_23_TH);
diff_check_6 = simplify(diff_check_6)


% D_q(3,1)
term26 = D_q(3,1)
yy(1,27) = cos(q2 + q4);
term =  subs(term26, yy(1,27), w);
TH27 = diff(term,w);
term27 = simplify(term-TH27*w)

yy(1,28) = sin(q2 + q4);
term =  subs(term27, yy(1,28), w);
TH28 = diff(term,w);
term28 = simplify(term-TH28*w)

D_31_TH = (TH27*yy(1,27)+TH28*yy(1,28));
diff_check_7 = eval(D_q(3,1)-D_31_TH);
diff_check_7 = simplify(diff_check_7)

% D_q(3,2)
term28 = D_q(3,2)
yy(1,29) = cos(q4);
term =  subs(term28, yy(1,29), w);
TH29 = diff(term,w);
term29 = simplify(term-TH29*w)

yy(1,30) = sin(q4);
term =  subs(term29, yy(1,30), w);
TH30 = diff(term,w);
term30 = simplify(term-TH30*w)

yy(1,31) = 1;
TH31 = m4*a^2 - 2*m4*a*xc4 + m4*xc4^2 + m4*zc4^2 + I4yy;
term31 = term30-TH31

D_32_TH = (TH29*yy(1,29)+TH30*yy(1,30)+TH31);
diff_check_8 = eval(D_q(3,2)-D_32_TH);
diff_check_8 = simplify(diff_check_8)

% D_q(3,3)
term31 = D_q(3,3)
yy(1,32) = 1;
TH32 = m4*a^2 - 2*m4*a*xc4 + m4*xc4^2 + m4*zc4^2 + I4yy;
term32 = term31-TH32

D_33_TH = TH32;
diff_check_9 = eval(D_q(3,3)-D_33_TH)

% gq(1,1)
gq(1,1)
G_11_TH = 0;

% gq(2,1)
term32 = gq(2,1)
yy(1,33) = cos(q2)*cos(q4);
term =  subs(term32, yy(1,33), w);
TH33 = diff(term,w);
term33 = simplify(term-TH33*w)

yy(1,34) = cos(q2)*sin(q4);
term =  subs(term33, yy(1,34), w);
TH34 = diff(term,w);
term34 = simplify(term-TH34*w)

yy(1,35) = cos(q4)*sin(q2);
term =  subs(term34, yy(1,35), w);
TH35 = diff(term,w);
term35 = simplify(term-TH35*w)

yy(1,36) = sin(q2)*sin(q4);
term =  subs(term35, yy(1,36), w);
TH36 = diff(term,w);
term36 = simplify(term-TH36*w)

yy(1,37) = cos(q2);
term =  subs(term36, yy(1,37), w);
TH37 = diff(term,w);
term37 = simplify(term-TH37*w)

yy(1,38) = sin(q2);
term =  subs(term37, yy(1,38), w);
TH38 = diff(term,w);
term38 = simplify(term-TH38*w)

G_21_TH = (TH33*yy(1,33)+TH34*yy(1,34)+TH35*yy(1,35)+TH36*yy(1,36)+TH37*yy(1,37)+TH38*yy(1,38));
diff_check_10 = eval(gq(2,1)-G_21_TH);
diff_check_10 = simplify(diff_check_10)

% gq(3,1)
term38 = gq(3,1)
yy(1,39) = cos(q2)*cos(q4);
term =  subs(term38, yy(1,39), w);
TH39 = diff(term,w);
term39 = simplify(term-TH39*w)

yy(1,40) = sin(q2)*sin(q4);
term =  subs(term39, yy(1,40), w);
TH40 = diff(term,w);
term40 = simplify(term-TH40*w)

yy(1,41) = sin(q2 + q4);
term =  subs(term40, yy(1,41), w);
TH41 = diff(term,w);
term41 = simplify(term-TH41*w)

G_31_TH = (TH39*yy(1,39)+TH40*yy(1,40)+TH41*yy(1,41));
diff_check_11 = eval(gq(3,1)-G_31_TH);
diff_check_11 = simplify(diff_check_11)


D_th = [D_11_TH, D_12_TH, D_13_TH ;
        D_21_TH, D_22_TH, D_23_TH ;
        D_31_TH, D_32_TH, D_33_TH ];
    
g_th = [G_11_TH; G_21_TH; G_31_TH ];
TH = [TH1; TH2; TH3; TH4; TH5; TH6; TH7; TH8; TH9; TH10; TH11; TH12; TH13; TH14; TH15; TH16; TH17; TH18; TH19; TH20; TH21; TH22; TH23; TH24; TH25; TH26; TH27; TH28; TH29; TH30; TH31; TH32; TH33; TH34; TH35; TH36; TH37; TH38; TH39; TH40; TH41];



%% 8,9. Finding unique parameters from M and G
TH_unique = unique(TH);
TH_unique_len = length(TH_unique)
% 


%% 10 Finding Coriolis matrix 
for i=1:3
    for j=1:3
        for k=1:3
            c_th(i,j,k)=(1/2)*((diff(D_th(k,j),q(i)))+(diff(D_th(k,i),q(j)))-(diff(D_th(i,j),q(k))));
        end
    end
end

for k=1:3
    for j=1:3
        C_th(k,j)=c_th(1,j,k)*q1dot+c_th(2,j,k)*q2dot+c_th(3,j,k)*q4dot;
    end
end



%% 11. Verify that the parameterization
D_q_TH_difference = simplify(D_q-D_th)
C_q_TH_difference = simplify(C-C_th)
g_q_TH_difference = simplify(gq-g_th)



%% 12 finding regressor
syms q1ddot q2ddot q4ddot
qddot = [q1ddot;q2ddot;q4ddot];

syms TH1 TH2 TH3 TH4 TH5 TH6 TH7 TH8 TH9 TH10 
syms TH11 TH12 TH13 TH14 TH15 TH16 TH17 TH18 TH19 TH20 
syms TH21 TH22 TH23 TH24 TH25 TH26 TH27 TH28 TH29 TH30 
syms TH31 TH32 TH33 TH34 TH35 TH36 TH37 TH38 TH39 TH40 
syms TH41

TH_syms = [TH1 TH2 TH3 TH4 TH5 TH6 TH7 TH8 TH9 TH10 TH11 TH12 TH13 TH14 TH15 TH16 TH17 TH18 TH19 TH20 TH21 TH22 TH23 TH24 TH25 TH26 TH27 TH28 TH29 TH30 TH31 TH32 TH33 TH34 TH35 TH36 TH37 TH38 TH39 TH40 TH41].';

D_11_TH = (TH1+TH2*yy(1,2)+TH3*yy(1,3)+TH4*yy(1,4)+TH5*yy(1,5)+TH6*yy(1,6)+TH7*yy(1,7)+TH8*yy(1,8)+TH9*yy(1,9)+TH10);
D_12_TH = (TH11*yy(1,11)+TH12*yy(1,12)+TH13*yy(1,13)+TH14*yy(1,14));
D_13_TH = (TH15*yy(1,15)+TH16*yy(1,16));
D_21_TH = (TH17*yy(1,17)+TH18*yy(1,18)+TH19*yy(1,19)+TH20*yy(1,20));
D_22_TH = (TH21*yy(1,21)+TH22*yy(1,22)+TH23);
D_23_TH = (TH24*yy(1,24)+TH25*yy(1,25)+TH26);
D_31_TH = (TH27*yy(1,27)+TH28*yy(1,28));
D_32_TH = (TH29*yy(1,29)+TH30*yy(1,30)+TH31);
D_33_TH = TH32;
G_11_TH = 0;
G_21_TH = (TH33*yy(1,33)+TH34*yy(1,34)+TH35*yy(1,35)+TH36*yy(1,36)+TH37*yy(1,37)+TH38*yy(1,38));
G_31_TH = (TH39*yy(1,39)+TH40*yy(1,40)+TH41*yy(1,41));

D_th = [D_11_TH, D_12_TH, D_13_TH ;
        D_21_TH, D_22_TH, D_23_TH ;
        D_31_TH, D_32_TH, D_33_TH ];
    
g_th = [G_11_TH; G_21_TH; G_31_TH ];


for i=1:3
    for j=1:3
        for k=1:3
            c_th(i,j,k)=(1/2)*((diff(D_th(k,j),q(i)))+(diff(D_th(k,i),q(j)))-(diff(D_th(i,j),q(k))));
        end
    end
end

for k=1:3
    for j=1:3
        C_th(k,j)=c_th(1,j,k)*q1dot+c_th(2,j,k)*q2dot+c_th(3,j,k)*q4dot;
    end
end


aux = D_th*qddot+C_th*qdot+g_th;

for i=1:3
    for j=1:41
        Y(i,j) = diff(aux(i),TH_syms(j));
    end
end


aux2 = D_q*qddot+C*qdot+gq;

% check if regreesor and original equation match
aux_YTHETA_difference = simplify(aux2-Y*TH)

%% 13. evaluating TH


%gravity
g=9.81;

%Local CM coordinates

xc1=-0.00443422;
yc1=0.12189039;
zc1=-0.00066489;
xc2=-0.00236983;
yc2=0.03105614;
zc2=0.01542114;
xc3=-0.03825858;
yc3=0.20750770;
zc3=0.00003309;
xc4=0.01095471;
yc4=-0.00002567;
zc4=0.14053900;

%Moments of inertia
I1xx=0.29486350;
I1xy=-0.00795023;
I1xz=-0.00009311;
I1yy=0.11350017;
I1yz=-0.00018711;
I1zz=0.25065343;

I2xx=0.02606840;
I2xy=-0.00001346;
I2xz=-0.00011701;
I2yy=0.01472202;
I2yz=0.00003659;
I2zz=0.01934814;

I3xx=0.13671601;
I3xy=-0.01680434;
I3xz=0.00000510;
I3yy=0.00588354;
I3yz=-0.00000530;
I3zz=0.13951371;

I4xx=0.03952350;
I4xy=0.00000189;
I4xz=0.00003117;
I4yy=0.04008214;
I4yz=0.00000131;
I4zz=0.00210299;

%Link masses
m1=10.76768767;
m2=3.87493756;
m3=1.80228141;
m4=1.06513649; %elbow+blank link

%Other kinematic parameters
a=0.045;
d=0.35;

TH = eval(TH)



sim('hw5_sim');

plot(time,q_1,'k')
hold on
plot(time,q_2,'r')
hold on
plot(time,q_3,'b')
hold on
grid on
xlabel('time')
ylabel('q_i')
legend('q1','q2','q3')

