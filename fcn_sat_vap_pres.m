function es = fcn_sat_vap_pres(T)
%% saturation vapor pressure
% a function of temperature in Kelvin
% given by Fanning & Weaver (1996)
% units are millibars
%%
es = 6.112.*exp(17.67*(T-273.15)./(T-29.65));

end

