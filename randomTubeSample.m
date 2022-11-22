function [M,POmega,Omega,A] = randomTubeSample( T, sizeOmega, Omega )

temp = size(T);
n = [temp(1),temp(2)];

if sizeOmega > prod(n)
    error('非法采样数目')
end

M = zeros(size(T));
POmega = zeros(n);
A = zeros(size(T));

if exist('Omega','var')==0    
    idx = randi( prod(n), sizeOmega, 1 );
    Omega = unique(idx);
end

while length(Omega) < sizeOmega
    idx = [ Omega; randi( prod(n) , sizeOmega-length(Omega), 1 )]; 
    Omega = unique(idx);
end
    
Omega = sort( Omega(1:sizeOmega) );

[i,j] = ind2sub( n, Omega );
subs = [i,j];

for k=1:length(Omega)
    M(subs(k,1),subs(k,2),:) = T(subs(k,1),subs(k,2),:);
    POmega(subs(k,1),subs(k,2)) = 1;
    A(subs(k,1),subs(k,2),:) = 1;
end


