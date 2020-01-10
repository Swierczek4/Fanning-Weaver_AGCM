function QLW = fcn_QLW(T,epsP,sigma)
%% infrared emission
% a function of Temperature in Kelvin and other parameters
% given by Fanning & Weaver (1996)
% units are W/m/m
%%
QLW = sigma.*epsP.*T.^4;

end

