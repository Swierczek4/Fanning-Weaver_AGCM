function rel = fcn_rel_humid(T,Q)
%% relative humidity
% a function of temperature and humidity
% given by Fanning & Weaver (1996)
% unitless (kg/kg)
%%

rel = Q./fcn_sat_spec_humid(T);

end

