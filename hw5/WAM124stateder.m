function zdot = WAM124stateder(z,tau,TH)

z1 = z(1:3); 
z2 = z(4:6);

q1 = z(1);
q2 = z(2);
q4 = z(3);

q1dot = z(4);
q2dot = z(5);
q4dot = z(6);

TH1 = TH(1); TH2 = TH(2); TH3 = TH(3); TH4 = TH(4); TH5 = TH(5); TH6 = TH(6); TH7 = TH(7); TH8 = TH(8); TH9 = TH(9); TH10 = TH(10);
TH11 = TH(11); TH12 = TH(12); TH13 = TH(13); TH14 = TH(14); TH15 = TH(15); TH16 = TH(16); TH17 = TH(17); TH18 = TH(18); TH19 = TH(19); TH20 = TH(20);
TH21 = TH(21); TH22 = TH(22); TH23 = TH(23); TH24 = TH(24); TH25 = TH(25); TH26 = TH(26); TH27 = TH(27); TH28 = TH(28); TH29 = TH(29); TH30 = TH(30);
TH31 = TH(31); TH32 = TH(32); TH33 = TH(33); TH34 = TH(34); TH35 = TH(35); TH36 = TH(36); TH37 = TH(37); TH38 = TH(38); TH39 = TH(39); TH40 = TH(40);
TH41 = TH(41); 

D_th = [ TH1 + TH10 + TH2*cos(2*q2 + 2*q4) + TH3*sin(2*q2 + 2*q4) + TH6*cos(q4) + TH9*sin(q4) + TH7*cos(2*q2 + q4) + TH8*sin(2*q2 + q4) + TH4*cos(2*q2) + TH5*sin(2*q2), TH11*cos(q2 + q4) + TH13*sin(q2 + q4) + TH12*cos(q2) + TH14*sin(q2), TH15*cos(q2 + q4) + TH16*sin(q2 + q4);
                                                                                            TH17*cos(q2 + q4) + TH19*sin(q2 + q4) + TH18*cos(q2) + TH20*sin(q2),                                  TH23 + TH21*cos(q4) + TH22*sin(q4),    TH26 + TH24*cos(q4) + TH25*sin(q4);
                                                                                                                          TH27*cos(q2 + q4) + TH28*sin(q2 + q4),                                  TH31 + TH29*cos(q4) + TH30*sin(q4),                                  TH32
        ];
g_th =[
                                                                                                                       0;
 TH37*cos(q2) + TH38*sin(q2) + TH33*cos(q2)*cos(q4) + TH34*cos(q2)*sin(q4) + TH35*cos(q4)*sin(q2) + TH36*sin(q2)*sin(q4);
                                                         TH41*sin(q2 + q4) + TH39*cos(q2)*cos(q4) + TH40*sin(q2)*sin(q4)
       ];

C_th = [
q4dot*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + (TH9*cos(q4))/2 - (TH6*sin(q4))/2 + (TH8*cos(2*q2 + q4))/2 - (TH7*sin(2*q2 + q4))/2) + q2dot*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + TH8*cos(2*q2 + q4) - TH7*sin(2*q2 + q4) + TH5*cos(2*q2) - TH4*sin(2*q2)), q2dot*(TH13*cos(q2 + q4) - TH11*sin(q2 + q4) + TH14*cos(q2) - TH12*sin(q2)) + q4dot*((TH13*cos(q2 + q4))/2 + (TH16*cos(q2 + q4))/2 - (TH11*sin(q2 + q4))/2 - (TH15*sin(q2 + q4))/2) + q1dot*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + TH8*cos(2*q2 + q4) - TH7*sin(2*q2 + q4) + TH5*cos(2*q2) - TH4*sin(2*q2)), q1dot*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + (TH9*cos(q4))/2 - (TH6*sin(q4))/2 + (TH8*cos(2*q2 + q4))/2 - (TH7*sin(2*q2 + q4))/2) + q2dot*((TH13*cos(q2 + q4))/2 + (TH16*cos(q2 + q4))/2 - (TH11*sin(q2 + q4))/2 - (TH15*sin(q2 + q4))/2) + q4dot*(TH16*cos(q2 + q4) - TH15*sin(q2 + q4));
                                     q4dot*((TH19*cos(q2 + q4))/2 - (TH28*cos(q2 + q4))/2 - (TH17*sin(q2 + q4))/2 + (TH27*sin(q2 + q4))/2) - q1dot*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + TH8*cos(2*q2 + q4) - TH7*sin(2*q2 + q4) + TH5*cos(2*q2) - TH4*sin(2*q2)),                                                                                     q4dot*((TH22*cos(q4))/2 - (TH21*sin(q4))/2) - q1dot*((TH13*cos(q2 + q4))/2 - (TH19*cos(q2 + q4))/2 - (TH11*sin(q2 + q4))/2 + (TH17*sin(q2 + q4))/2 + (TH14*cos(q2))/2 - (TH20*cos(q2))/2 - (TH12*sin(q2))/2 + (TH18*sin(q2))/2),                                                                                                         q2dot*((TH22*cos(q4))/2 - (TH21*sin(q4))/2) + q4dot*(TH25*cos(q4) - TH24*sin(q4)) - q1dot*((TH16*cos(q2 + q4))/2 - (TH19*cos(q2 + q4))/2 - (TH15*sin(q2 + q4))/2 + (TH17*sin(q2 + q4))/2);
                       - q1dot*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + (TH9*cos(q4))/2 - (TH6*sin(q4))/2 + (TH8*cos(2*q2 + q4))/2 - (TH7*sin(2*q2 + q4))/2) - q2dot*((TH19*cos(q2 + q4))/2 - (TH28*cos(q2 + q4))/2 - (TH17*sin(q2 + q4))/2 + (TH27*sin(q2 + q4))/2),                                                                                                                                                               - q2dot*((TH22*cos(q4))/2 - (TH21*sin(q4))/2) - q1dot*((TH13*cos(q2 + q4))/2 - (TH28*cos(q2 + q4))/2 - (TH11*sin(q2 + q4))/2 + (TH27*sin(q2 + q4))/2),                                                                                                       - q2dot*((TH25*cos(q4))/2 - (TH30*cos(q4))/2 - (TH24*sin(q4))/2 + (TH29*sin(q4))/2) - q1dot*((TH16*cos(q2 + q4))/2 - (TH28*cos(q2 + q4))/2 - (TH15*sin(q2 + q4))/2 + (TH27*sin(q2 + q4))/2)
                       
                       ];


D_th_eval = (D_th);
g_th_eval = (g_th);
C_th_eval = (C_th);

z2dot = D_th_eval\(-C_th_eval*z2-g_th_eval+tau);

zdot = [z2;z2dot];







% zdot = [.1;.1;.1;.1;.1;.1];


% yy_eval(1,1) = 1;
% yy_eval(1,2) = cos(2*q2 + 2*q4);
% yy_eval(1,3) = sin(2*q2 + 2*q4);
% yy_eval(1,4) = cos(2*q2);
% yy_eval(1,5) = sin(2*q2);
% yy_eval(1,6) = cos(q4);
% yy_eval(1,7) = cos(2*q2 + q4);
% yy_eval(1,8) = sin(2*q2 + q4);
% yy_eval(1,9) = sin(q4);
% yy_eval(1,10) = 1;
% yy_eval(1,11) = cos(q2 + q4);
% yy_eval(1,12) = cos(q2);
% yy_eval(1,13) = sin(q2 + q4);
% yy_eval(1,14) = sin(q2);
% yy_eval(1,15) = cos(q2 + q4);
% yy_eval(1,16) = sin(q2 + q4);
% yy_eval(1,17) = cos(q2 + q4);
% yy_eval(1,18) = cos(q2);
% yy_eval(1,19) = sin(q2 + q4);
% yy_eval(1,20) = sin(q2);
% yy_eval = eval(yy);

% D_11_TH = (TH1+TH2*yy_eval(1,2)+TH3*yy_eval(1,3)+TH4*yy_eval(1,4)+TH5*yy_eval(1,5)+TH6*yy_eval(1,6)+TH7*yy_eval(1,7)+TH8*yy_eval(1,8)+TH9*yy_eval(1,9)+TH10);
% D_12_TH = (TH11*yy_eval(1,11)+TH12*yy_eval(1,12)+TH13*yy_eval(1,13)+TH14*yy_eval(1,14));
% D_13_TH = (TH15*yy_eval(1,15)+TH16*yy_eval(1,16));
% D_21_TH = (TH17*yy_eval(1,17)+TH18*yy_eval(1,18)+TH19*yy_eval(1,19)+TH20*yy_eval(1,20));
% D_22_TH = (TH21*yy_eval(1,21)+TH22*yy_eval(1,22)+TH23);
% D_23_TH = (TH24*yy_eval(1,24)+TH25*yy_eval(1,25)+TH26);
% D_31_TH = (TH27*yy_eval(1,27)+TH28*yy_eval(1,28));
% D_32_TH = (TH29*yy_eval(1,29)+TH30*yy_eval(1,30)+TH31);
% D_33_TH = TH32;
