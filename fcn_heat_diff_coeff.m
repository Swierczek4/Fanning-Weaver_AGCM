function nu = fcn_heat_diff_coeff(lat)
%% heat diffusion coefficient
% a function of latitude in degrees
% given by Fanning & Weaver (1996)
% units are m*m/s
%%
rad = lat.*pi./180;
nu = 3e6.*(0.81-1.08*(sin(rad)).^2 + ...
    0.74.*(sin(rad)).^4);

end

