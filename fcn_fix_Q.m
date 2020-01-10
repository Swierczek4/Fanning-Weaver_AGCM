function Q = fcn_fix_Q(T,Q)
%% humidity correction
% a function of Temperature in Kelvin and other parameters
% given by Fanning & Weaver (1996)
% units are m/yr
%%
rel = fcn_rel_humid(T,Q);
rel=rel>=.85;
qs = fcn_sat_spec_humid(T);

[n,m] = size(rel);
for ii=1:n
    for jj=1:m
        if (rel(ii,jj)==1)
           Q(ii,jj) = 0.85*qs(ii,jj); 
        end
    end
end

end


