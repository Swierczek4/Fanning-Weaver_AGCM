function EVAP = fcn_evap(T,Q,rho_air,rho_sea,SST,U,qsatSST)
%% evaporation
% a function of Temperature in Kelvin and other parameters
% given by Fanning & Weaver (1996)
% units are m/yr
%%
EVAP = (rho_air/rho_sea)*31536000.*fcn_Ce(T,SST,U).*U.*(qsatSST-Q);

end

