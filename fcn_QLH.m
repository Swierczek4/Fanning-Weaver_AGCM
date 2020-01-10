function QLH = fcn_QLH(T,Lnu,rho_sea,rho_air,Hq,Q,dt)
%% sensible heat flux for atmosphere
% a function of pressure
% given by Fanning & Weaver (1996)
% units are W/m/m
%%
QLH = rho_sea*Lnu.*fcn_precip(T,rho_air,rho_sea,Hq,Q,dt)./31536000;

end

