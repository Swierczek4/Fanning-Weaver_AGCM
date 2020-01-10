function epsP = fcn_plan_emissivity(lat)
%% atmospheric emissivity
% a function of latitude in degrees
% given by Fanning & Weaver (1996)
% unitless
%%
sv = sin(lat.*pi./180);

epsP = 0.5531 - 0.1296.*sv + 0.6796.*sv.^2 + ...
    0.7116.*sv.^3 - 2.7940.*sv.^4 - 1.3592.*sv.^5 + ...
    3.8831.*sv.^6 + 0.8348.*sv.^7 - 1.9536.*sv.^8;

end

