function [PRECIP,Q] = fcn_precip(T,rho_air,rho_sea,Hq,Q,dt)
%% precipitation
% a function of Temperature in Kelvin and other parameters
% given by Fanning & Weaver (1996)
% units are m/yr
%%
rel = fcn_rel_humid(T,Q);
qs = fcn_sat_spec_humid(T);

PRECIP = rho_air*Hq*31536000.*fcn_H(rel).*(Q -...
    .85.*qs)./rho_sea./dt;

[n,m] = size(rel);
for ii=1:n
    for jj=1:m
        if (rel(ii,jj)==1)
           Q(ii,jj) = 0.85*qs(ii,jj); 
        end
    end
end
end

