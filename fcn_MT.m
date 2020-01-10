function MT = fcn_MT(Q,rho_air,Hq,kappa,DX,DY,order)
%% eddy diffusive horizontal moisture transport parameterization
% a pde
% given by Fanning & Weaver (1996)
% units are W/m/m
%%
coeff = rho_air*Hq;
Qxx = pde_finite_diff_p(Q,1,2,order,DX);
Qyy = pde_finite_diff_n(Q,2,2,order,DY);
term1 = kappa.*(Qxx+Qyy);
kappay = pde_finite_diff_n(kappa,2,1,order,DY);
Qy = pde_finite_diff_n(Q,2,1,order,DY);
term2 = kappay.*Qy;
MT = coeff*(term1+term2);

% coeff = rho_air*Hq;
% Qx = pde_finite_diff_p(Q,1,1,order,DX);
% Qy = pde_finite_diff_n(Q,2,1,order,DY);
% term1 = kappa.*Qx;
% term2 = kappa.*Qy;
% term3 = pde_finite_diff_p(term1,1,1,order,DX);
% term4 = pde_finite_diff_n(term2,2,1,order,DY);
% MT = coeff*(term3+term4);

end

