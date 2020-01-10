function alpha = fcn_coalbedo(lat)
%% coalbedo
% a function of latitude in degrees
% given by Fanning & Weaver (1996)
% unitless
%%
alpha = 0.7995 - 0.315.*(sin(lat.*pi./180)).^2;
end

