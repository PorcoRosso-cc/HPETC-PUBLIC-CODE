function [Qu,Qv,Up,Up0,Vp,Vp0,lambda,beta] = produceWeightedTensor(T,n1,n2,n3,r,gaussianNoise)

[U, ~, V, ~] = tsvd(T,'full');

U = U(:,1:r,:);
V = V(:,1:r,:);

trank = r; 

if exist('gaussianNoise','var')==0
    gaussianNoise = 1e-2;
end

T0 = fft(T,[],3);
U0 = fft(U,[],3);
V0 = fft(V,[],3);

Up = zeros(n1,trank,n3);
Vp = zeros(n2,trank,n3);
Tg = zeros(n1,n2,n3);
Up0 = zeros(n1,n1-trank,n3);
Vp0 = zeros(n2,n2-trank,n3);
for i=1:n3
   %G = sqrt(gaussianNoise)*randn(size(n1,n2));
   G = sqrt(gaussianNoise)*randn(n1,n2);
   Tg(:,:,i) = T0(:,:,i) + G;
   [Ug, ~, Vg] = svd(Tg(:,:,i));
   Up(:,:,i) = Ug(:,1:trank);
   Vp(:,:,i) = Vg(:,1:trank);
   [qu,~] = qr(Up(:,:,i));
   Up0(:,:,i) = qu(:,trank+1:n1);
   [qv,~] = qr(Vp(:,:,i));
   Vp0(:,:,i) = qv(:,trank+1:n2);
end
%为了设置信任值参数，我们计算先验U V与真实U V之间的最大主角,单独编写一个M文件
Cu = zeros(1,n3);
Cv = zeros(1,n3);
for i=1:n3
    Cu(i) = cosofmatrix(U0(:,:,i),Up(:,:,i));
    Cv(i) = cosofmatrix(V0(:,:,i),Vp(:,:,i));
end
%根据最大主角，取tan()获得信任值参数
lambda = zeros(1,n3);
beta = zeros(1,n3);
for i=1:n3
   lambda(i) = min(sqrt(tan(Cu(i))),1);
   beta(i) = min(sqrt(tan(Cv(i))),1);
end
%构造权重矩阵Qu Qv
Qu = zeros(n1,n1,n3);
Qv = zeros(n2,n2,n3);

for i=1:n3
    Qu(:,:,i) = lambda(i)*(Up(:,:,i)*Up(:,:,i)')+Up0(:,:,i)*Up0(:,:,i)';
    Qv(:,:,i) = beta(i)*(Vp(:,:,i)*Vp(:,:,i)')+Vp0(:,:,i)*Vp0(:,:,i)';
end
    


% Qu = ifft(Qu,[],3);
% Qv = ifft(Qv,[],3);


