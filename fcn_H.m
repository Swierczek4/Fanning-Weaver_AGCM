function H = fcn_H(rel)
%% relative humidity switch
% a function of relative humidity
% given by Fanning & Weaver (1996)
% unitless
%%
H = rel>=0.85;
H = double(H);
end

