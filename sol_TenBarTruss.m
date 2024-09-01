function [Q, stress, R, K] = sol_TenBarTruss(r)
  
    % 定義各參數數值  
    
    %節點座標
    node = [18.28 9.14;18.28 0;9.14 9.14;9.14 0;0 9.14;0 0];
    ec=[3 5;1 3;4 6;2 4;3 4;1 2;4 5;3 6;2 3;1 4];
    ed=[5 6 9 10; 1 2 5 6; 7 8 11 12; 3 4 7 8; 5 6 7 8; 1 2 3 4; 7 8 9 10; 5 6 11 12; 3 4 5 6;  1 2 7 8]; 
    E = 200 *10^9;
    A(1:6) = pi*r(1)^2;
    A(7:10) = pi*r(2)^2;
    L(1:6) = 9.14;
    L(7:10) = 9.14*sqrt(2);     
    % 開一個空白的剛性矩陣 (stiffness matrix)
    K = zeros(12);          

    % 計算 stiffness matrix     
    for i=1:10
        C=(node(ec(i,2),1)-node(ec(i,1),1))/L(i); 
        S=(node(ec(i,2),2)-node(ec(i,1),2))/L(i);
        k=(A(i)*E/L(i)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;-C*C -C*S C*C C*S; -C*S -S*S C*S S*S]);        
        ev=ed(i,:);
        for x=1:4
           for y=1:4
                K( ev(1,x), ev(1,y) ) = K( ev(1,x), ev(1,y) ) + k(x,y);
           end
        end
    end   

    % 建立力矩陣   
    F(4) = -10^7 ;
    F(8) = -10^7 ;

    % 計算位移量 (F = KQ)
    Q = inv(K(1:8,1:8))*F';
    Q(12)=0;      
    

    % 計算應力 (stress)
    for i=1:10
        C=(node(ec(i,2),1)-node(ec(i,1),1))/L(i); 
        S=(node(ec(i,2),2)-node(ec(i,1),2))/L(i);
        stress(i)=(E/L(i))*[-C -S C S]*Q((ed(i,:)));  
    end  

    % (optional) compute reactions
    R = K(9:12,1:12) * Q;
end