


function [x, P] = mu_g(x, P, yacc, Ra, g0)

hx = Qq(x)'*g0;
[Q0, Q1, Q2, Q3] = dQqdq(x);
J0 = Q0'*g0;
J1 = Q1'*g0;
J2 = Q2'*g0;
J3 = Q3'*g0;
Hx = [J0, J1, J2, J3];
S = Hx*P*Hx' + Ra;
K = P*Hx'/S;
x = x + K*(yacc - hx);
P = P - K*S*K';

end