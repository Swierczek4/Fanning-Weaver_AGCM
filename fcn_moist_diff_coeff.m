function kappa = fcn_moist_diff_coeff(lat)
%% moisture diffusion coefficient
% a function of latitude in degrees
% given by Fanning & Weaver (1996)
% units are m*m/s
%%
sv = abs(sin(lat.*pi./180));

kappa = 1.7e6.*(1.9823 - 17.3501.*sv + 117.2489.*sv.^2 - ...
    274.1129.*sv.^3 + 258.2244.*sv.^4 - 85.7967.*sv.^5);

end

