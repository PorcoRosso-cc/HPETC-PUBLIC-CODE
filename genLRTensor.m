function T = genLRTensor(n1, n2, n3, r)
% 生成低秩张量

% n1 = 50; n2 = 50; n3 = 100; %张量维数
% r = 5; %张量秩

%第一步，先生成随机矩阵
u1 = randn(n1, r); u2 = randn(n2, r); u3 = randn(n3, r);
% [U1 temp] = qr(u1); [U2 temp] = qr(u2); [U3 temp] = qr(u3);
% U1 = U1(:,1:r); U2 = U2(:,1:r); U3 = U3(:,1:r);
U1 = u1(:,1:r); U2 = u2(:,1:r); U3 = u3(:,1:r);
for i=1:n3
    T(:,:,i) = U1*diag(U3(i,:))*U2';
end
