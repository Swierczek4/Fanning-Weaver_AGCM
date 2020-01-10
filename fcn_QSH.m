function QSH = fcn_QSH(SST,T,rho_air,Ch,Crhoa,U)
%% sensible heat flux
% a function of temperature in Kelvin, other parameters
% and U, scalar surface wind speed
% given by Fanning & Weaver (1996)
% units are W/m/m
%%
QSH = rho_air*Crhoa*Ch.*U.*(SST-T);

end

