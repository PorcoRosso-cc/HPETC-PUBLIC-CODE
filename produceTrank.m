 function [trank] = produceTrank(s1)
%PRODUCETRANK 此处显示有关此函数的摘要
%   此处显示详细说明

len = length(s1);
tol = 0;

for i=1:len
    tol = tol + s1(i);
end

trank = 1;
num = 0;
while true
    num = num + s1(trank);
    if num/tol > 0.6
        break;
    end
    trank = trank + 1;
end

end

