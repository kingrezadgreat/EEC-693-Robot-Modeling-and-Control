function u=robust_passivity(t,z,zd)

g=9.8;
q1=z(1);
q2=z(2);
q4=z(3);
q1dot=z(4);
q2dot=z(5);
q4dot=z(6);

q1d=zd(1);
q2d=zd(2);
q4d=zd(3);
q1ddot=zd(4);
q2ddot=zd(5);
q4ddot=zd(6);
q1dddot=zd(7);
q2dddot=zd(8);
q4dddot=zd(9);

qtilde=[q1-q1d;q2-q2d;q4-q4d];
qtildedot=[q1dot-q1ddot;q2dot-q2ddot;q4dot-q4ddot];

%Gains
L=0.0001*diag([1 2 1]);
K=0.0001*diag([1 2 1]);

v=[q1ddot;q2ddot;q4ddot]-L*qtilde;
a=[q1dddot;q2dddot;q4dddot]-L*qtildedot;
r=qtildedot+L*qtilde;

a1 = a(1);
a2 = a(2);
a3 = a(3);

v1 = v(1);
v2 = v(2);
v3 = v(3);

Yav = [a1, a1*cos(2*q2 + 2*q4) - v1*(q2d*sin(2*q2 + 2*q4) + q4d*sin(2*q2 + 2*q4)) - q1d*v2*sin(2*q2 + 2*q4) - q1d*v3*sin(2*q2 + 2*q4), a1*sin(2*q2 + 2*q4) + v1*(q2d*cos(2*q2 + 2*q4) + q4d*cos(2*q2 + 2*q4)) + q1d*v2*cos(2*q2 + 2*q4) + q1d*v3*cos(2*q2 + 2*q4), a1*cos(2*q2) - q1d*v2*sin(2*q2) - q2d*v1*sin(2*q2), a1*sin(2*q2) + q1d*v2*cos(2*q2) + q2d*v1*cos(2*q2), a1*cos(2*q2 + q4) - v1*(q2d*sin(2*q2 + q4) + (q4d*sin(2*q2 + q4))/2) - q1d*v2*sin(2*q2 + q4) - (q1d*v3*sin(2*q2 + q4))/2, a1*sin(2*q2 + q4) + v1*(q2d*cos(2*q2 + q4) + (q4d*cos(2*q2 + q4))/2) + q1d*v2*cos(2*q2 + q4) + (q1d*v3*cos(2*q2 + q4))/2, a1*cos(q4) - (q1d*v3*sin(q4))/2 - (q4d*v1*sin(q4))/2, a1*sin(q4) + (q1d*v3*cos(q4))/2 + (q4d*v1*cos(q4))/2, a2*cos(q2 + q4) - v3*(q2d*sin(q2 + q4) + q4d*sin(q2 + q4)) - v2*(q2d*sin(q2 + q4) + q4d*sin(q2 + q4)) + a3*cos(q2 + q4), v2*(q2d*cos(q2 + q4) + q4d*cos(q2 + q4)) + v3*(q2d*cos(q2 + q4) + q4d*cos(q2 + q4)) + a2*sin(q2 + q4) + a3*sin(q2 + q4), a2*cos(q2) - q2d*v2*sin(q2), a2*sin(q2) + q2d*v2*cos(q2),  0,                                                    0,                                                    0,       0,            0,            0,       0,       0, q1d, sign(q1d),   0,   0;
       0,                                                                                                    q1d*v1*sin(2*q2 + 2*q4),                                                                                                   -q1d*v1*cos(2*q2 + 2*q4),                                   q1d*v1*sin(2*q2),                                  -q1d*v1*cos(2*q2),                                                                                                    q1d*v1*sin(2*q2 + q4),                                                                                                   -q1d*v1*cos(2*q2 + q4),                          a3*cos(q4) - q4d*v3*sin(q4),                          a3*sin(q4) + q4d*v3*cos(q4),                                                                                                         a1*cos(q2 + q4),                                                                                                         a1*sin(q2 + q4),                  a1*cos(q2),                  a1*sin(q2), a2, a2*sin(q4) + (q2d*v3*cos(q4))/2 + (q4d*v2*cos(q4))/2, a2*cos(q4) - (q2d*v3*sin(q4))/2 - (q4d*v2*sin(q4))/2,      a3, cos(q2 + q4), sin(q2 + q4), cos(q2), sin(q2),   0,         0, q2d,   0;
       0,                                                                                                    q1d*v1*sin(2*q2 + 2*q4),                                                                                                   -q1d*v1*cos(2*q2 + 2*q4),                                                  0,                                                  0,                                                                                                (q1d*v1*sin(2*q2 + q4))/2,                                                                                               -(q1d*v1*cos(2*q2 + q4))/2,                      a2*cos(q4) + (q1d*v1*sin(q4))/2,                      a2*sin(q4) - (q1d*v1*cos(q4))/2,                                                                                                         a1*cos(q2 + q4),                                                                                                         a1*sin(q2 + q4),                           0,                           0,  0,                                  -(q2d*v2*cos(q4))/2,                                   (q2d*v2*sin(q4))/2, a2 + a3, cos(q2 + q4), sin(q2 + q4),       0,       0,   0,         0,   0, q4d];
 
TH_0 =[
    0.4693
   -0.0099
   -0.0048
   -0.1137
    0.0158
   -0.0540
   -0.0060
    0.0508
    0.0194
    0.0000
   -0.0000
   -0.0037
   -0.0006
    0.3612
    0.0389
    0.1015
    0.1179
   -0.3557
    1.4685
    0.4993
    6.7627
    1.4385
    2.0414
    4.5858
    1.3665
]*.1;

TH_max = TH_0*.1;

deadzone=0.1;
rho=8;
% rho=2.78;
% rho=5.6158;
Y = Yav;
s=Y'*r;
% size(s)
% size(TH_0)
if norm(s)>deadzone
    dTh=-rho*s/norm(s);
else
    dTh=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
end
Th_hat=TH_0+dTh;

u=Y*Th_hat-K*r;







% %In-class example: 2-link planar RR manipulator
% g=9.8;
% q1=z(1);
% q2=z(2);
% q1dot=z(3);
% q2dot=z(4);
% 
% q1d=zd(1);
% q2d=zd(2);
% q1ddot=zd(3);
% q2ddot=zd(4);
% q1dddot=zd(5);
% q2dddot=zd(6);
% 
% qtilde=[q1-q1d;q2-q2d];
% qtildedot=[q1dot-q1ddot;q2dot-q2ddot];
% 
% %Gains
% L=1*diag([1 1]);
% K=1*diag([1 1]);
% 
% v=[q1ddot;q2ddot]-L*qtilde;
% a=[q1dddot;q2dddot]-L*qtildedot;
% r=qtildedot+L*qtilde;
% 
% %Form the regressor
% Y(1,1)=a(1);
% Y(1,2)=cos(q2)*(2*a(1)+a(2))-sin(q2)*q2dot*v(1)-sin(q2)*(q1dot+q2dot)*v(2);
% Y(1,3)=a(2);
% Y(1,4)=g*cos(q1);
% Y(1,5)=g*cos(q1+q2);
% Y(2,1)=0;
% Y(2,2)=cos(q2)*a(1)+sin(q2)*q1dot*v(1);
% Y(2,3)=a(1)+a(2);
% Y(2,4)=0;
% Y(2,5)=Y(1,5);
% 
% % TH_0=[42.1;0.96;12.32;7.8;3.2];
% TH_0=[4.1;5.96;1.32;.8;.2];
% 
% deadzone=0.1;
% rho=32.78;
% % rho=2.78;
% % rho=5.6158;
% s=Y'*r;
% if norm(s)>deadzone
%     dTh=-rho*s/norm(s);
% else
%     dTh=[0;0;0;0;0];
% end
% Th_hat=TH_0+dTh;
% 
% u=Y*Th_hat-K*r;
