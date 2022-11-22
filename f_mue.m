function [out] = f_mue(mue, L, R, Upre, Upre_orth, Vpre, Vpre_orth, l, b)
%F_MUE 此处显示有关此函数的摘要
%   此处显示详细说明
out = 1/2 * trace( ( l^2/(mue*l^2+1)^2 * Upre + 1/(mue+1)^2 * Upre_orth  )*(L*L')...
        + ( b^2/(mue*b^2+1)^2 * Vpre + 1/(mue+1)^2 * Vpre_orth )*(R*R') );
end

