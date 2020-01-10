function S = fcn_shortwave_dist(lat)
%% annual distribution of shortwave radiation
% a function of latitude in degrees
% given by Fanning & Weaver (1996)
% unitless
%%
S = 1.5.*(1-(sin(lat.*pi./180).^2));

end

