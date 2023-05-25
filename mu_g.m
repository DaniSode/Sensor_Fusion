function [x, P] = mu_g(x, P, yacc, Ra, g0)

[hx,Hx] = h(x);
S = Hx*P*Hx' + R;
K = P*Hx'/S;
x = x + K*(y - hx);
P = (eye(size(K,1)) - K*Hx) * P;