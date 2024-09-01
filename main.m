clc;
clear;
obj = @object;
nonlcon = @nonlcon;
sol_TenBarTruss = @sol_TenBarTruss;
%----------GA algorithm Demo----------%
r0  =  [0.1,0.1];
A = [];% 輸入線性不等式拘束條件的係數矩陣A與B
 % A*X<B
b = [];
Aeq = []; % 輸入線性等式的拘束條件係數矩陣Aeq與Beq
beq = []; % Aeq*X=Beq
LB = [0.001;0.001]; % 設計變數的下限
UB = [0.5;0.5]; % 設計變數的上限
options = optimoptions('fmincon','Display','final','Algorithm','sqp'); % 演算法的參數使用內設值
[r,fval,exitflag] = fmincon(obj,r0,A,b,Aeq,beq,LB,UB,nonlcon,options);
[Q, stress, R, K] = sol_TenBarTruss(r);
% x:最佳解
