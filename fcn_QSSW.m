function QSSW = fcn_QSSW(So,S,Coal,Co)
%% incoming shortwave radiation source
% a function of solar constant and other parameters
% given by Fanning & Weaver (1996)
% units are W/m/m
%%
QSSW = So.*S.*Coal.*(1-Co)./4;

end

