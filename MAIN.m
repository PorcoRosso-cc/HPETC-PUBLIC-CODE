clear;
clc;

% load datasets: PlanetLab
load_path = 'data/PlanetLab490';
load(load_path);
Tensor = PlanetLab;
[n1,n2,n3] = size(Tensor);
Tensor = Tensor ./ max(max(max(Tensor)));
R = 13;

fprintf('%15s %13s %13s %13s %13s %13s \n', 'Algorithm','SampleRatio' ,'RSE', 'MAE', 'RMSE', 'Time');
Omega = [];
for i =1:19
    p = 0.05 * i; 
    % Algorithm:HPETC 1e-6
    t = tic;
    [T_hpetc, Omega1, POmega] = PWTNN(Tensor, p, R, 1e-6, Omega);
    time = toc(t);
    % Save result
    rse = norm(T_hpetc(:)-Tensor(:))/norm(Tensor(:));
    mae = sum(abs(T_hpetc - Tensor), 'all')/(n1*n2*n3);
    rmse = norm(T_hpetc(:)-Tensor(:))/sqrt(n1*n2*n3);
    fprintf('%15s %13.2f %13e %13e %13e %13d \n','HPETC Gau=1e-6', ...
        p, rse, mae, rmse, time);
end

