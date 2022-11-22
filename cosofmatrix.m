function c = cosofmatrix(A,B)
%求解最大主角
[~,n2] = size(A);
c = 0;
ca = zeros(1,n2);
for i=1:n2
   minangle = 100;
   minj = 1;
   [~,m] = size(B);
   for j=1:m
      temp = subspace(A(:,i),B(:,j));
      if temp<minangle
         minangle = temp;
         minj = j;
      end
   end
   ca(i) = minangle;
   B(:,minj) = [];
end
c = max(ca);