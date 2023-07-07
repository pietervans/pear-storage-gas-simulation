function c =  getInitialC(Ku, Kv, fu, fv, tr, CONSTS, varargin)
% Return initial c for Newton iterations by using linear approximation of R_u and R_v
% Inputs:
%  - Ku, Kv: matrices K_u and K_v from assignment
%            (from integrals 1 and 3 in Equations (5) and (6))
%  - fu, fv: vectors f_u and f_v from assignment
%            (from integrals 1 and 3 in Equations (5) and (6))
%  - tr: triangulation object
%  - CONSTS: struct
%  - varargin: optional boolean describing whether or not to compute Cv

% Add contributions from integral 2 in assignment (Equations (5) and (6))
[Kv_temp, fv_temp] = get_Kvfv2(tr, CONSTS);

c1 = (Ku+get_Ku2(tr, CONSTS))\fu;

if length(varargin) == 1
    compute_cv = varargin{1};
else
    compute_cv = true;
end

if compute_cv
    c2 = Kv\(fv+fv_temp-Kv_temp*c1);
else
    c2 = zeros(size(c1));
end
c = [c1; c2];
end
