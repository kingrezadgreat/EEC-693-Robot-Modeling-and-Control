function J = cost_function(X,L1,L2,t)

% a = X(1:12);
% b = X(13:24);
%
% a = reshape(a,[4,3]);
% b = reshape(b,[4,3]);

a1 = X(1:4);
a2 = X(5:8);
a4 = X(9:12);

b1 = X(13:16);
b2 = X(17:20);
b4 = X(21:24);

w = 2*pi/10;
J = 0;
Y_svd = 0;
Y_cond = 0;
for tt = t
    
    q1=0;q2=0;q4=0;
    q1d=0;q2d=0;q4d=0;
    q1dd=0;q2dd=0;q4dd=0;
    
    %         q_temp = 0;
    %         qd_temp = 0;
    %         qdd_temp = 0;
    
    for j = 1:4
        q1 = q1 + a1(j)*sin(w*j*tt)/(w*j) - b1(j)*cos(w*j*tt)/(w*j);
        q2 = q2 + a2(j)*sin(w*j*tt)/(w*j) - b2(j)*cos(w*j*tt)/(w*j);
        q4 = q4 + a4(j)*sin(w*j*tt)/(w*j) - b4(j)*cos(w*j*tt)/(w*j);
        
        q1d = q1d + a1(j)*cos(w*j*tt) + b1(j)*sin(w*j*tt);
        q2d = q2d + a2(j)*cos(w*j*tt) + b2(j)*sin(w*j*tt);
        q4d = q4d + a4(j)*cos(w*j*tt) + b4(j)*sin(w*j*tt);
        
        q1dd = q1dd - a1(j)*w*j*sin(w*j*tt) - b1(j)*w*j*cos(w*j*tt);
        q2dd = q2dd - a2(j)*w*j*sin(w*j*tt) - b2(j)*w*j*cos(w*j*tt);
        q4dd = q4dd - a4(j)*w*j*sin(w*j*tt) - b4(j)*w*j*cos(w*j*tt);
        
        %             q_temp   = q_temp +             a(i,j)*sin(w*(i)*tt)/(w*(i)) -           b(i,j)*cos(w*(i)*tt)/(w*(i)) ;
        %             qd_temp  = qd_temp  +   a(i,j)*(w*(i))*cos(w*(i)*tt)/(w*(i)) +   b(i,j)*(w*(i))*sin(w*(i)*tt)/(w*(i)) ;
        %             qdd_temp = qdd_temp - a(i,j)*(w*(i))^2*sin(w*(i)*tt)/(w*(i)) + b(i,j)*(w*(i))^2*cos(w*(i)*tt)/(w*(i)) ;
    end
    
    q = [q1;q2;q4];
    qd = [q1d;q2d;q4d];
    qdd = [q1dd;q2dd;q4dd];
    %         q(j,1) = q_temp;
    %         qd(j,1) = qd_temp;
    %         qdd(j,1) = qdd_temp;
    
    Y = regressor(q,qd,qdd);
    Y_svd = Y_svd + min(svd(Y));
    Y_cond = Y_cond + cond(Y);
    
end
J = L1*(Y_cond./100) + L2/(Y_svd/100);

end