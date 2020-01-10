function Ch = fcn_Ch(T,SST,U)
%% Stanton number
% a function of temperature, SST, and wind speed
% given by Fanning & Weaver (1996)
% unitless
%%
Ch = 0.94.*fcn_Ce(T,SST,U);

end

