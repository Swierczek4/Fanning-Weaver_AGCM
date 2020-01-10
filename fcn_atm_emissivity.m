function epsA = fcn_atm_emissivity(lat)
%% atmospheric emissivity
% a function of latitude in degrees
% given by Fanning & Weaver (1996)
% unitless
%%
sv = sin(lat.*pi./180);

epsA = 0.8666 + 0.0408.*sv - 0.2553.*sv.^2 - ...
    0.4660.*sv.^3 + 0.9877.*sv.^4 + 2.0257.*sv.^5 - ...
    2.3374.*sv.^6 - 3.1990.*sv.^7 + 2.8581.*sv.^8 + ...
    1.6070.*sv.^9 - 1.2685.*sv.^10;

end

