function [L,R,r_norm] = mySPG(POmega, M_fi, tau, Up, Uo, Vp, Vo, l, b, L, R)
% 求解LASSO子问题
% 设置的参数：1. step length: alpha_min, alpha_max, aplpha,
%            2. descent parameter: gama:(0,1)
%            3. 初始化变量: X = Project[X], r = M_fi - POmega.*X, g = r
%            4. integer linesearch history length M >=1, （先不设置）


alpha_max = 1;
alpha = 1; % alpha初始值


X = L * R';
r = POmega .* X - M_fi;
g = r;

N = 0; %迭代轮数

theta = norm(r,'fro');
eps = 1e-4;

while true   
    while true
       L_temp = L - alpha * g * R;
       R_temp = R - alpha * g' * L;
       
       [L_new, R_new] = projectedToFeasibleSet(tau, L_temp, R_temp, Up, Uo, Vp, Vo, l, b);
       X_new = L_new * R_new';
       r_new = POmega .* X_new - M_fi;
       
       if norm(r_new,'fro') < norm(r,'fro')
           break;
       else
           alpha = alpha * 0.5;
       end
       if alpha < 1e-5
           break;
       end
    end
    theta_old = theta;
    theta = norm(r_new,'fro');
%==========experiment_three======    
    if abs(theta-theta_old)/theta_old < eps
        break;
    end
%=============================
    L = L_new; 
    R = R_new; 
    r = r_new; 
    g = r;
    alpha = alpha_max;
    N = N + 1;
    if N>10
        break;
    end
end
r_norm = norm(r,'fro');
end

