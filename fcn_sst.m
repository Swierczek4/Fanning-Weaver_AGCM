function SST = fcn_sst(sst)
%% convert the sst file into Kelvin except for -999 values
logic = sst>-990;
SST = sst + 273.15;
SST = logic.*SST;
end

