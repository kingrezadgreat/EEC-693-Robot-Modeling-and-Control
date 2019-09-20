function Yav = regressor(q,qd,qdd)

q1 = q(1);
q2 = q(2);
q4 = q(3);

q1d = qd(1);
q2d = qd(2);
q4d = qd(3);

q1dd = qdd(1);
q2dd = qdd(2);
q4dd = qdd(3);

Y = [q1dd, q1dd*cos(2*q2 + 2*q4) - 2*q1d*q2d*sin(2*q2 + 2*q4) - 2*q1d*q4d*sin(2*q2 + 2*q4), q1dd*sin(2*q2 + 2*q4) + 2*q1d*q2d*cos(2*q2 + 2*q4) + 2*q1d*q4d*cos(2*q2 + 2*q4), q1dd*cos(2*q2) - 2*q1d*q2d*sin(2*q2), q1dd*sin(2*q2) + 2*q1d*q2d*cos(2*q2), q1dd*cos(2*q2 + q4) - 2*q1d*q2d*sin(2*q2 + q4) - q1d*q4d*sin(2*q2 + q4), q1dd*sin(2*q2 + q4) + 2*q1d*q2d*cos(2*q2 + q4) + q1d*q4d*cos(2*q2 + q4),   q1dd*cos(q4) - q1d*q4d*sin(q4),   q1dd*sin(q4) + q1d*q4d*cos(q4), - sin(q2 + q4)*q2d^2 - 2*sin(q2 + q4)*q2d*q4d - sin(q2 + q4)*q4d^2 + q2dd*cos(q2 + q4) + q4dd*cos(q2 + q4), cos(q2 + q4)*q2d^2 + 2*cos(q2 + q4)*q2d*q4d + cos(q2 + q4)*q4d^2 + q2dd*sin(q2 + q4) + q4dd*sin(q2 + q4), - sin(q2)*q2d^2 + q2dd*cos(q2), cos(q2)*q2d^2 + q2dd*sin(q2),    0,                              0,                              0,           0,            0,            0,       0,       0;
    0,                                                          q1d^2*sin(2*q2 + 2*q4),                                                         -q1d^2*cos(2*q2 + 2*q4),                      q1d^2*sin(2*q2),                     -q1d^2*cos(2*q2),                                                    q1d^2*sin(2*q2 + q4),                                                   -q1d^2*cos(2*q2 + q4),   - sin(q4)*q4d^2 + q4dd*cos(q4),     cos(q4)*q4d^2 + q4dd*sin(q4),                                                                                          q1dd*cos(q2 + q4),                                                                                        q1dd*sin(q2 + q4),                   q1dd*cos(q2),                 q1dd*sin(q2), q2dd, q2dd*sin(q4) + q2d*q4d*cos(q4), q2dd*cos(q4) - q2d*q4d*sin(q4),        q4dd, cos(q2 + q4), sin(q2 + q4), cos(q2), sin(q2);
    0,                                                          q1d^2*sin(2*q2 + 2*q4),                                                         -q1d^2*cos(2*q2 + 2*q4),                                    0,                                    0,                                                (q1d^2*sin(2*q2 + q4))/2,                                               -(q1d^2*cos(2*q2 + q4))/2, (sin(q4)*q1d^2)/2 + q2dd*cos(q4), q2dd*sin(q4) - (q1d^2*cos(q4))/2,                                                                                          q1dd*cos(q2 + q4),                                                                                        q1dd*sin(q2 + q4),                              0,                            0,    0,             -(q2d^2*cos(q4))/2,              (q2d^2*sin(q4))/2, q2dd + q4dd, cos(q2 + q4), sin(q2 + q4),       0,       0];

Yf = [q1d, sign(q1d), 0 ,0;
    0, 0, q2d, 0;
    0, 0, 0, q4d];

Yav = [Y,Yf];

end
