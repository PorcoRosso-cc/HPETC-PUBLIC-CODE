function [out] = grads_f_mue(mue, L, R, Upre, Upre_orth, Vpre, Vpre_orth, l, b)
%GRADS_F_MUE 此处显示有关此函数的摘要
%   此处显示详细说明
out = trace( ( (-2*l^4)/(mue*l^2+1)^3 * Upre + (-2)/(mue+1)^3 * Upre_orth )*(L*L')... 
        + ( (-2*b^4)/(mue*b^2+1)^3 * Vpre + (-2)/(mue+1)^3 * Vpre_orth )*(R*R') );
end

