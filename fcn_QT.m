function QT = fcn_QT(T,rho_air,Ha,Crhoa,nu,DX,DY,order)
%% eddy diffusive horizontal heat transport parameterization
% a pde
% given by Fanning & Weaver (1996)
% units are W/m/m
%%
coeff = rho_air*Ha*Crhoa;
Txx = pde_finite_diff_p(T,1,2,order,DX);
Tyy = pde_finite_diff_n(T,2,2,order,DY);
term1 = nu.*(Txx+Tyy);
nuy = pde_finite_diff_n(nu,2,1,order,DY);
Ty = pde_finite_diff_n(T,2,1,order,DY);
term2 = nuy.*Ty;
QT = coeff*(term1+term2);

% coeff = rho_air*Ha*Crhoa;
% Tx = pde_finite_diff_p(T,1,1,order,DX);
% Ty = pde_finite_diff_n(T,2,1,order,DY);
% term1 = nu.*Tx;
% term2 = nu.*Ty;
% term3 = pde_finite_diff_p(term1,1,1,order,DX);
% term4 = pde_finite_diff_n(term2,2,1,order,DY);
% QT = coeff*(term3+term4);

end

