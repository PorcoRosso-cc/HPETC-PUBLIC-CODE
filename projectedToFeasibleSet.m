function [L_new, R_new] = projectedToFeasibleSet(tau, L, R, Up, Uo, Vp, Vo, l, b)
% ��XͶӰ�����н⼯�����X_fea
% ��ͨ����� f(mue)<=tau, ����� mue ֵ�� Ȼ���� L, R
mue = 1;
fmue = f_mue(mue, L, R, Up, Uo, Vp, Vo, l, b);
grads_fmue = grads_f_mue(mue, L, R, Up, Uo, Vp, Vo, l, b);
converged = false;
N = 0;
maxIter = 100;
% if fmue < tau
%    L_new = L;
%    R_new = R;
%    return
% end

while ~converged
   N = N + 1;
   mue_new = mue - (fmue-tau)/grads_fmue;
   fmue_new = f_mue(mue_new, L, R, Up, Uo, Vp, Vo, l, b);
   
   mue_old = mue;
   mue = mue_new;
   fmue = fmue_new;
   grads_fmue = grads_f_mue(mue, L, R, Up, Uo, Vp, Vo, l, b);
   
%===========�˹�����=========
%    if abs(mue_old-mue) < 1e-6
%        break;
%    end
%    if fmue_new  < tau
%        break;
%    end
%=============================

%===========��ʵ����=========
   if abs(mue_old-mue)/abs(mue_old) < 1e-2
       break;
   end
   if fmue_new  < tau
       break;
   end
    if N>maxIter
        break;
    end
%=============================

    
end

L_new = ( 1/(mue*l^2+1) * Up + 1/(mue+1) * Uo ) * L;
R_new = ( 1/(mue*b^2+1) * Vp + 1/(mue+1) * Vo ) * R;

end

