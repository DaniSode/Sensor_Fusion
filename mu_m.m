


function [x, P] = mu_m(x, P, mag, m0, Rm)

hx = Qq(x)'*m0;
[Q0, Q1, Q2, Q3] = dQqdq(x);
J0 = Q0'*m0;
J1 = Q1'*m0;
J2 = Q2'*m0;
J3 = Q3'*m0;
Hx = [J0, J1, J2, J3];
S = Hx*P*Hx' + Rm;
K = P*Hx'/S;
x = x + K*(mag - hx);
P = P - K*S*K';

end