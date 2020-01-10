function Ce = fcn_Ce(T,SST,U)
%% Dalton number
% a function of temperature, SST, and wind speed
% given by Fanning & Weaver (1996)
% unitless
%%
Ce = 1e-3.*(1.0022-0.822.*(T-SST)+.0266.*U);
[nn,mm] = size(Ce);

for ii=1:nn
    for jj=1:mm
        if Ce(ii,jj)<6e-5
            Ce(ii,jj) = 6e-5;
        elseif Ce(ii,jj)>2.19e-3
            Ce(ii,jj) = 2.19e-3;
        end
    end
end

end

