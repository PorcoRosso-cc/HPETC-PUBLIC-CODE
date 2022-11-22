function [U,S,V]=t_SVD(Tensor)
%================== 张量tSVD分解 ===================%
%输入：Tensor 数据张量
%输出：U V 特征空间 S对角
%==================================================%
[I,J,K]=size(Tensor);
T=fft(Tensor,[],3); %傅里叶变换
U=zeros(I,I,K);  S=zeros(I,J,K); V=zeros(J,J,K);
for k=1:K
    [u,s,v]=svd(T(:,:,k)); %对第k个面做svd分解
    U(:,:,k)=u;
    S(:,:,k)=s;
    V(:,:,k)=v;
end
U=ifft(U,[],3); % 这里分块得到相应的子矩阵 A1,A2
S=ifft(S,[],3); % 
V=ifft(V,[],3);