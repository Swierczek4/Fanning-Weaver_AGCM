%% Make a grid of solar scattering coefficient

load('file_mask.mat')

[n,m] = size(mask);
Co = 0.65*ones(n,m);
land = zeros(n,m);

for ii=1:n
    for jj=1:m
        if mask(ii,jj)==0
            Co(ii,jj) = 0.3;
            land(ii,jj) = 1;
        end
    end
end

sea = mask;

clear n m ii jj mask
