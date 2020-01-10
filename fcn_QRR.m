function QRR = fcn_QRR(SST,T,epsO,epsA,sigma)
%% radiative flux
% a function of temperature in Kelvin
% given by Fanning & Weaver (1996)
% units are W/m/m
%%
QRR = sigma*(epsO.*SST.^4 - epsA.*T.^4);

end

