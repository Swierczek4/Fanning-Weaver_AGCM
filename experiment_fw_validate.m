clear
close all
clc

tic()

%% preliminaries & fixed parameters
acc_settings        % file has custom colors, colormaps, and plot settings
dt = 3600/4;
num_steps = floor(17280000/dt);
dx = 5;
dy = 4;
start_lat = -90;
stop_lat = 90;
start_lon = -180;
stop_lon = 175;
% parameters
fcn_solar_scattering_coeff
radius = 6371000;
rho_air = 1.25;
rho_sea = 1024;
epsO = 0.96;
Ha = 8400;
Hq = 1800;
Lnu = 2.5e6;
So = 1360;
sigma = 5.67e-8;
Crhoa = 1000;
%% end preliminaries & fixed parameters

%% grid
LON_GRID = start_lon:dx:stop_lon;
LAT_GRID = start_lat:dy:stop_lat;
[LON_MESH,LAT_MESH] = meshgrid(LON_GRID,LAT_GRID);
Nx = length(LON_GRID);
Ny = length(LAT_GRID);
DY = dy*pi*radius/180;
DX = dx*pi*radius.*abs(cos(pi/180.*LAT_MESH))./180;
LAT_MESH = flipud(LAT_MESH);
%% end grid

%% read netcdf and matlab files
LONTP = ncread('file_precip.nc','LONT');
LATTP = ncread('file_precip.nc','LATT');
PRECIP_TRUE = ncread('file_precip.nc','PRECIP');
LONTS = ncread('file_sat.nc','LONT');
LATTS = ncread('file_sat.nc','LATT');
SAT_TRUE = ncread('file_sat.nc','SAT');
load('file_sst.mat');
load('file_wind.mat');
SST = fcn_sst(sst);
U = wind;
clear sst wind
load('output_Q.mat');
load('output_T.mat');
load('output_PRECIP_Time_Series');
load('output_T_Time_Series');
load('output_mean_PRECIP.mat');

%% end read netcdf and matlab files

%% fixed location-based parameters
Coal = fcn_coalbedo(LAT_MESH);
epsA = fcn_atm_emissivity(LAT_MESH);
nu = fcn_heat_diff_coeff(LAT_MESH);
kappa = fcn_moist_diff_coeff(LAT_MESH);
epsP = fcn_plan_emissivity(LAT_MESH);
S = fcn_shortwave_dist(LAT_MESH);
qsatSST = fcn_sat_spec_humid(SST);
%%

%% temperature plot
cm = acc_colormap('es_coolwarm');

figure()
set(gcf, 'Position', [1, 1, 1400, 700])
colormap(cm)
contourf(LON_GRID,LAT_GRID,output_T-273.15,'LineStyle','none','LevelList',z);
colorbar('eastoutside');
caxis(coloraxis)
hold on
contour(LON_GRID,LAT_GRID,0+land,'Color','k')
axis(coords)
xtickformat('degrees')
ytickformat('degrees')
title('surface air temperature estimate [deg C] after 200 days',...
    'FontWeight','Normal','FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
acc_plots
hold off
print('-r200','plot_SAT_EST','-dpng')
%% end temperature plot

%% true temperature plot
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
colormap(cm)
contourf(LON_GRID,LAT_GRID,SAT_TRUE','LineStyle','none','LevelList',z4);
colorbar('eastoutside');
caxis(coloraxis4)
hold on
contour(LON_GRID,LAT_GRID,0+land,'Color','k')
axis(coords)
xtickformat('degrees')
ytickformat('degrees')
title('climatological surface air temperature [deg C]',...
    'FontWeight','Normal','FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
acc_plots
hold off
print('-r200','plot_SAT_TRUE','-dpng')
%% end true temperature plot

%% error plot for temperature
cm = acc_colormap('cmo_balance');

figure()
set(gcf, 'Position', [1, 1, 1400, 700])
colormap(cm)
contourf(LON_GRID,LAT_GRID,output_T-273.15-SAT_TRUE','LineStyle','none','LevelList',z1);
colorbar('eastoutside');
caxis(coloraxis1)
hold on
contour(LON_GRID,LAT_GRID,0+land,'Color','k')
axis(coords)
xtickformat('degrees')
ytickformat('degrees')
title('surface air temperature error [deg C] (estimate - climatology)',...
    'FontWeight','Normal','FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
acc_plots
hold off
print('-r200','plot_SAT_ERR','-dpng')
%% end error plot for temperature

%% precipitation plot
cm = acc_colormap('cmo_rain');
PRECIP_EST = fcn_precip(output_T,rho_air,rho_sea,Hq,output_Q,dt);
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
colormap(cm)
contourf(LON_GRID,LAT_GRID,output_mean_PRECIP,'LineStyle','none','LevelList',z3);
colorbar('eastoutside');
caxis(coloraxis3)
hold on
contour(LON_GRID,LAT_GRID,0+land,'Color','k')
axis(coords)
xtickformat('degrees')
ytickformat('degrees')
title('precipitation estimate [m/yr] after 200 days',...
    'FontWeight','Normal','FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
acc_plots
hold off
print('-r200','plot_PRECIP_EST','-dpng')
%% end precipitation plot

%% precipitation plot
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
colormap(cm)
contourf(LON_GRID,LAT_GRID,0.365.*PRECIP_TRUE','LineStyle','none','LevelList',z6);
colorbar('eastoutside');
caxis(coloraxis6)
hold on
contour(LON_GRID,LAT_GRID,0+land,'Color','k')
axis(coords)
xtickformat('degrees')
ytickformat('degrees')
title('climatological precipitation [m/yr]','FontWeight','Normal','FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
acc_plots
hold off
print('-r200','plot_PRECIP_TRUE','-dpng')
%% end true precipitation plot

%% error plot for precipitation
cm = acc_colormap('cmo_balance');
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
colormap(cm)
contourf(LON_GRID,LAT_GRID,output_mean_PRECIP-0.365.*PRECIP_TRUE','LineStyle','none','LevelList',z2);
colorbar('eastoutside');
caxis(coloraxis2)
hold on
contour(LON_GRID,LAT_GRID,0+land,'Color','k')
axis(coords)
xtickformat('degrees')
ytickformat('degrees')
title('precipitation error [m/yr] (estimate - climatology)',...
    'FontWeight','Normal','FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
acc_plots
hold off
print('-r200','plot_PRECIP_ERR','-dpng')
%% end error plot for precipitation

%% plot global mean surface temp time series
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
plot((1:num_steps)*dt/86400,output_T_Time_Series,'Color',Color(:,color1),'LineWidth',lw)
title('global mean surface temperature',...
    'FontWeight','Normal','FontSize',16)
xlabel('time (days)')
ylabel('degrees C')
acc_plots
print('-r200','plot_Time_Series_T','-dpng')
%% end plot global mean surface temp time series

%% plot global mean precipitation time series
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
plot((1:num_steps)*dt/86400,output_PRECIP_Time_Series,'Color',Color(:,color2),'LineWidth',lw)
title('global mean precipitation',...
    'FontWeight','Normal','FontSize',16)
xlabel('time (days)')
ylabel('m/yr')
acc_plots
print('-r200','plot_Time_Series_PRECIP','-dpng')
%% end plot global mean precipitation time series

%% plot of zonal mean temperature comparison
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
h1 = plot(LAT_GRID,mean(output_T-273.15,2),'Color',Color(:,color1),'LineWidth',lw);
hold on
h2 = plot(LAT_GRID,mean(SAT_TRUE,1),'Color',Color(:,color2),'LineWidth',lw);
xlim([-90 90])
title('zonal mean surface air temperature [deg C]',...
    'FontWeight','Normal','FontSize',16)
xlabel('Latitude')
ylabel('degrees C')
legend([h1(1),h2(1)],'estimate','climatological')
acc_plots
hold off
print('-r200','plot_SAT_ZONAL','-dpng')
%% end plot of zonal mean temperature comparison

%% plot of zonal mean precipitation comparison
figure()
set(gcf, 'Position', [1, 1, 1400, 700])
h1 = plot(LAT_GRID,mean(output_mean_PRECIP,2),'Color',Color(:,color1),'LineWidth',lw);
hold on
h2 = plot(LAT_GRID,mean(0.365*PRECIP_TRUE,1),'Color',Color(:,color2),'LineWidth',lw);
xlim([-90 90])
title('zonal mean precipitation [m/yr]','FontWeight','Normal','FontSize',16)
xlabel('Latitude')
ylabel('m/yr')
legend([h1(1),h2(1)],'estimate','climatological')
acc_plots
hold off
print('-r200','plot_PRECIP_ZONAL','-dpng')
%% end plot of zonal mean precipitation comparison

toc()