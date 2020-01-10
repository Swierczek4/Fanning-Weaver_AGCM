function qs = fcn_sat_spec_humid(T)
%% saturation specific humidity
% a function of temperature in Kelvin
% given by Fanning & Weaver (1996)
% unitless (kg/kg)
%%
es = fcn_sat_vap_pres(T);
qs = 0.622.*es./(1013.26 - 0.378.*es);

end

