function zdot=stateder(t,z,u)


n=length(z);
z_1=z(1:n/2);
z_2=z(n/2+1:n);

q1 = z(1);
q2 = z(2);
q4 = z(3);

q1d = z(4);
q2d = z(5);
q4d = z(6);

cm11 = [-0.00443422;0.12189039;-0.00066489;1];
cm22 = [-0.00236983;0.03105614;0.01542114;1];
cm33 = [-0.03825858;0.20750770;0.00003309;1];
cm44 = [0.01095471;-0.00002567;0.14053900;1];

g = 9.81;

% Local CM coordinates

xc1 = -0.00443422;
yc1 = 0.12189039;
zc1 = -0.00066489;
xc2 = -0.00236983;
yc2 = 0.03105614;
zc2 = 0.01542114;
xc3 = -0.03825858;
yc3 = 0.20750770;
zc3 = 0.00003309;
xc4 = 0.01095471;
yc4 = -0.00002567;
zc4 = 0.14053900;

%Moments of inertia taken at the CG, aligned with local frame
I1xx=0.13488033; I1xy = -0.00213041; I1xz = -0.00012485;
I1yx = -0.00213041; I1yy = 0.11328369; I1yz = 0.00068555;
I1zx = -0.00012485; I1zy = 0.00068555; I1zz = 0.09046330;

I2xx = 0.02140958 ;I2xy = 0.00027172 ;I2xz = 0.00002461;
I2yx = 0.00027172 ;I2yy = 0.01377875 ;I2yz = -0.00181920;
I2zx = 0.00002461 ;I2zy = -0.00181920 ;I2zz = 0.01558906;


I3xx = 0.05911077 ;I3xy = -0.00249612 ;I3xz = 0.00000738;
I3yx = -0.00249612 ;I3yy = 0.00324550 ;I3yz = -0.00001767;
I3zx = 0.00000738 ;I3zy = -0.00001767 ;I3zz = 0.05927043;

I4xx = 0.00093220 ;I4xy = 0.00000258 ;I4xz = -0.00031408;
I4yx = 0.00000258 ;I4yy = 0.00120070 ;I4yz = 0.00000227;
I4zx = -0.00031408 ;I4zy = 0.00000227 ;I4zz = 0.00084605;

% Link masses

m1 = 10.76768767;
m2 = 3.87493756;
m3 = 1.80228141;
m4 = 1.06513649; % elbow+blank link

% Other kinematic parameters

a = 0.045;
d = 0.35;

%Reflected drive inertias in kg-m^2
Jdrive1=0.205190;
Jdrive2=0.094428;
Jdrive4=0.094428;

TH = [I2xx/2 + I3xx/2 + I4xx/2 + I1yy + I3yy/2 + I2zz/2 + I4zz/2 + (a^2*m3)/2 + a^2*m4 + (d^2*m3)/2 + (d^2*m4)/2 + m1*xc1^2 + (m2*xc2^2)/2 + (m3*xc3^2)/2 + (m4*xc4^2)/2 + m2*yc2^2 + (m3*yc3^2)/2 + m4*yc4^2 + m1*zc1^2 + (m2*zc2^2)/2 + m3*zc3^2 + (m4*zc4^2)/2 + a*m3*xc3 - a*m4*xc4 - d*m3*yc3;
                                                                                                                                                                                                                        (m4*a^2)/2 - m4*a*xc4 + (m4*xc4^2)/2 - I4xx/2 + I4zz/2 - (m4*zc4^2)/2;
                                                                                                                                                                                                                                                                 m4*xc4*zc4 - a*m4*zc4 - I4xz;
                                                                                                                      I3yy/2 - I3xx/2 - I2xx/2 + I2zz/2 + (a^2*m3)/2 + (a^2*m4)/2 - (d^2*m3)/2 - (d^2*m4)/2 + (m2*xc2^2)/2 + (m3*xc3^2)/2 - (m3*yc3^2)/2 - (m2*zc2^2)/2 + a*m3*xc3 + d*m3*yc3;
                                                                                                                                                                                                                I3xy - I2xz + a*d*m3 + a*d*m4 - a*m3*yc3 + d*m3*xc3 - m3*xc3*yc3 + m2*xc2*zc2;
                                                                                                                                                                                                                                                                    -m4*(a^2 - xc4*a + d*zc4);
                                                                                                                                                                                                                                                                     m4*(a*zc4 - a*d + d*xc4);
                                                                                                                                                                                                                                                                   m4*(- a^2 + xc4*a + d*zc4);
                                                                                                                                                                                                                                                                     m4*(a*d + a*zc4 - d*xc4);
                                                                                                                                                                                                                                                                            I4yz - m4*yc4*zc4;
                                                                                                                                                                                                                                                                 m4*xc4*yc4 - a*m4*yc4 - I4xy;
                                                                                                                                                                                                                                  I2yz - I3yz - d*m4*yc4 - d*m3*zc3 - m2*yc2*zc2 + m3*yc3*zc3;
                                                                                                                                                                                                                                  a*m4*yc4 - I3xz - I2xy + a*m3*zc3 + m2*xc2*yc2 + m3*xc3*zc3;
                                                                                                                            I2yy + I4yy + I3zz + a^2*m3 + 2*a^2*m4 + d^2*m3 + d^2*m4 + m2*xc2^2 + m3*xc3^2 + m4*xc4^2 + m3*yc3^2 + m2*zc2^2 + m4*zc4^2 + 2*a*m3*xc3 - 2*a*m4*xc4 - 2*d*m3*yc3;
                                                                                                                                                                                                                                                           2*a*d*m4 + 2*a*m4*zc4 - 2*d*m4*xc4;
                                                                                                                                                                                                                                                                 2*m4*(- a^2 + xc4*a + d*zc4);
                                                                                                                                                                                                                                             m4*a^2 - 2*m4*a*xc4 + m4*xc4^2 + m4*zc4^2 + I4yy;
                                                                                                                                                                                                                                                                              -g*m4*(a - xc4);
                                                                                                                                                                                                                                                                                     g*m4*zc4;
                                                                                                                                                                                                                                                            g*(a*m3 + a*m4 + m2*xc2 + m3*xc3);
                                                                                                                                                                                                                                                            g*(d*m3 + d*m4 - m3*yc3 + m2*zc2)];
%%%End of TH listing  %%%%%
TH(1)=TH(1)+Jdrive1;
TH(14)=TH(14)+Jdrive2;
TH(17)=TH(17)+Jdrive4;

TH=[TH;1.4385;2.0414;4.5858;1.3665];
TH1 = TH(1);
TH2 = TH(2);
TH3 = TH(3);
TH4 = TH(4);
TH5 = TH(5);
TH6 = TH(6);
TH7 = TH(7);
TH8 = TH(8);
TH9 = TH(9);
TH10 = TH(10);
TH11 = TH(11);
TH12 = TH(12);
TH13 = TH(13);
TH14 = TH(14);
TH15 = TH(15);
TH16 = TH(16);
TH17 = TH(17);
TH18 = TH(18);
TH19 = TH(19);
TH20 = TH(20);
TH21 = TH(21);
TH22 = TH(22);
TH23 = TH(23);
TH24 = TH(24);
TH25 = TH(25);


M = [ TH1 + TH2*cos(2*q2 + 2*q4) + TH3*sin(2*q2 + 2*q4) + TH8*cos(q4) + TH9*sin(q4) + TH6*cos(2*q2 + q4) + TH7*sin(2*q2 + q4) + TH4*cos(2*q2) + TH5*sin(2*q2), TH10*cos(q2 + q4) + TH11*sin(q2 + q4) + TH12*cos(q2) + TH13*sin(q2), TH10*cos(q2 + q4) + TH11*sin(q2 + q4);
                                                                                     TH10*cos(q2 + q4) + TH11*sin(q2 + q4) + TH12*cos(q2) + TH13*sin(q2),                                  TH14 + TH16*cos(q4) + TH15*sin(q4),      TH17 + TH8*cos(q4) + TH9*sin(q4);
                                                                                                                   TH10*cos(q2 + q4) + TH11*sin(q2 + q4),                                    TH17 + TH8*cos(q4) + TH9*sin(q4),                                  TH17];
                                                                                                               
C = [q4d*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + (TH9*cos(q4))/2 - (TH8*sin(q4))/2 + (TH7*cos(2*q2 + q4))/2 - (TH6*sin(2*q2 + q4))/2) + q2d*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + TH7*cos(2*q2 + q4) - TH6*sin(2*q2 + q4) + TH5*cos(2*q2) - TH4*sin(2*q2)), q2d*(TH11*cos(q2 + q4) - TH10*sin(q2 + q4) + TH13*cos(q2) - TH12*sin(q2)) + q1d*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + TH7*cos(2*q2 + q4) - TH6*sin(2*q2 + q4) + TH5*cos(2*q2) - TH4*sin(2*q2)) + q4d*(TH11*cos(q2 + q4) - TH10*sin(q2 + q4)), q1d*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + (TH9*cos(q4))/2 - (TH8*sin(q4))/2 + (TH7*cos(2*q2 + q4))/2 - (TH6*sin(2*q2 + q4))/2) + q2d*(TH11*cos(q2 + q4) - TH10*sin(q2 + q4)) + q4d*(TH11*cos(q2 + q4) - TH10*sin(q2 + q4));
                                                                                                                                          -q1d*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + TH7*cos(2*q2 + q4) - TH6*sin(2*q2 + q4) + TH5*cos(2*q2) - TH4*sin(2*q2)),                                                                                                                                                                                                             q4d*((TH15*cos(q4))/2 - (TH16*sin(q4))/2),                                                                                                                                                         q4d*(TH9*cos(q4) - TH8*sin(q4)) + q2d*((TH15*cos(q4))/2 - (TH16*sin(q4))/2);
                                                                                                                              -q1d*(TH3*cos(2*q2 + 2*q4) - TH2*sin(2*q2 + 2*q4) + (TH9*cos(q4))/2 - (TH8*sin(q4))/2 + (TH7*cos(2*q2 + q4))/2 - (TH6*sin(2*q2 + q4))/2),                                                                                                                                                                                                            -q2d*((TH15*cos(q4))/2 - (TH16*sin(q4))/2),                                                                                                                                                                                                                                   0];
G = [ 0;
 TH18*cos(q2 + q4) + TH19*sin(q2 + q4) + TH20*cos(q2) + TH21*sin(q2);
                               TH18*cos(q2 + q4) + TH19*sin(q2 + q4)];                                                                                                                   



zdot=[q1d;q2d;q4d;inv(M)*(u-C*z_2-G)];





% 
% 
% %Two-link planar manipulator state derivatives
% 
% m1=10.9803;
% m2=8.5778;
% l1=0.6357;
% lc1=l1/2;  %to center of mass
% l2=0.4241;
% lc2=l2/2;  %to center of mass
% I1=25.406;
% I2=12.5467;
% g=9.8;
% 
% n=length(z);
% z_1=z(1:n/2);
% z_2=z(n/2+1:n);
% %find numerical values for matrices M, C and calculate state derivatives
% D(1,1)=m1*lc1^2+m2*(l1^2+lc2^2+2*l1*lc2*cos(z_1(2)))+I1+I2;
% D(1,2)=m2*(lc2^2+l1*lc2*cos(z_1(2)))+I2;
% D(2,1)=D(1,2);
% D(2,2)=m2*lc2^2+I2;
% 
% h=-m2*l1*lc2*sin(z_1(2));
% 
% C(1,1)=h*z_2(2);
% C(1,2)=h*z_2(2)+h*z_2(1);
% C(2,1)=-h*z_2(1);
% C(2,2)=0;
% 
% gg(1,1)=(m1*lc1+m2*l1)*g*cos(z_1(1))+m2*lc2*g*cos(z_1(1)+z_1(2));
% gg(2,1)=m2*lc2*g*cos(z_1(1)+z_1(2));
% 
% zdot=[z_2;inv(D)*(u-C*z_2-gg)];
% 
% 



