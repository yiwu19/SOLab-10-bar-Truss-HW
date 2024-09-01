function [g,geq]=nonlcon(r)
[Q, stress ] = sol_TenBarTruss(r);
g(1) = (Q(3,1)^2 + Q(4,1)^2)^0.5-0.02;
g(2) = max(abs(stress))-250*10^6;
geq = [];
end