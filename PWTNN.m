function [Tout,Omega, POmega] = PWTNN(Tensor,p, rank, Gau_var, Omega, Tensor_prior, lambda, beta)

% 1.construct weighted tensor W and Q 
[n1,n2,n3] = size(Tensor);
if exist('Gau_var','var') == 0
    Gau_var = 1e-6;
end
if exist('Omega','var') == 0
    Omega = [];
end
sizeOmega = ceil(n1*n2*p);
[M,POmega,Omega,~] = randomTubeSample(Tensor, sizeOmega,Omega);

if exist('Tensor_prior','var') == 0
    [Qu,Qv,Up,Up0,Vp,Vp0,lambda,beta] = produceWeightedTensor(Tensor,n1,n2,n3,rank,Gau_var);
else
    [Qu,Qv,Up,Up0,Vp,Vp0,lambda,beta] = produceWeightedTensorAPP(Tensor,Tensor_prior,n1,n2,n3,lambda,beta,rank);
end

% HPETC Al.
T_f = fft(Tensor,[],3); % target tensor, used to calculate error
M_f = fft(M,[],3); 
Tout_f = zeros(n1,n2,n3); 
iteration = 0;
for i=1:n3 
    tol = 1;
    if i>1
        tol = tol*1e2;
    end
    % rename paras
    T_fi = T_f(:,:,i);
    M_fi = M_f(:,:,i);
    Q = Qu(:,:,i); W = Qv(:,:,i);
    Upre = Up(:,:,i); Upre_orth = Up0(:,:,i);
    Vpre = Vp(:,:,i); Vpre_orth = Vp0(:,:,i);
    l = lambda(i); b = beta(i);
    N = 0;    
    % Initialize tau
    tau_true = norm_nuc(Q*T_fi*W);
    tau = tau_true * 1.1;
    % rename paras
    Upr = Upre*Upre'; Uo = Upre_orth*Upre_orth';
    Vpr = Vpre*Vpre'; Vo = Vpre_orth*Vpre_orth';   
    [U,S,V] = svd(M_fi,'econ');
    L = U * (S.^(1/2));
    R = V * (S.^(1/2));
    [L_old, R_old] = projectedToFeasibleSet(tau, L, R, Upr, Uo, Vpr, Vo, l, b);
    X = L_old * R_old';
    r_norm_old = norm(POmega .* X - M_fi, 'fro');
    roh = 0;
    gama = 1; % descent step
    while true
        iteration = iteration +1;
        N = N + 1;
        [L,R,r_norm] = mySPG(POmega, M_fi, tau, Upr, Uo, Vpr, Vo, l, b, L_old, R_old);
        if r_norm >= r_norm_old
            r_norm = r_norm_old;
            L = L_old;
            R = R_old;
            roh = roh + 1;
        else
            roh = 0;
        end
        if roh > 3
            break;
        end
        Tout_fi = L*R';
        if r_norm < tol
            break;
        end
        L_old = L; R_old = R; r_norm_old = r_norm;
        grads_LASSO = norm( pinv(Q)*(M_fi - POmega.*Tout_fi)*pinv(W), inf);
        tau_new = tau - gama*(r_norm-tol)/grads_LASSO;
        tau = tau_new;
        if tau < tau_true
            break;
        end
        if tau < 0
            break;
        end
    end
    Tout_f(:,:,i) = Tout_fi;
end
Tout = ifft(Tout_f,[],3);

end


