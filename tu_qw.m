function [x, P] = tu_qw(x, P, omega, T, Rw)

F = eye(size(x, 1)) + 1/2*Somega(omega)*T;
G = 1/2*Sq(x)*T;

x = F*x;
P = F*P*F' + G*Rw*G';

end
