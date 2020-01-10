clear
close all
clc

tic()

%% preliminaries & fixed parameters
acc_settings        % file has custom colors, colormaps, and plot settings
dt = 3600/7;
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
%% end grid

%% read netcdf and matlab files
load('file_sst.mat');
load('file_wind.mat');
SST = fcn_sst(sst);
U = wind;
clear sst wind
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

%% ICs
T_old = 273.15*ones(Ny,Nx);
output_T = T_old;
Q_old = zeros(Ny,Nx);
output_Q = Q_old;
output_T_Time_Series = zeros(num_steps,1);
output_PRECIP_Time_Series = zeros(num_steps,1);
[PRECIP,output_Q] = fcn_precip(output_T,rho_air,rho_sea,Hq,output_Q,dt);
output_mean_PRECIP = PRECIP;
%% end ICs

%% initial plot for movie
cm = acc_colormap('cmo_rain');

figure()
set(gcf, 'Position', [1, 1, 1600, 900])
colormap(cm)
contourf(LON_GRID,LAT_GRID,PRECIP,'LineStyle','none','LevelList',z3);
colorbar('eastoutside');
hold on
contour(LON_GRID,LAT_GRID,0+land,'Color','k')
caxis(coloraxis3)
axis(coords)
xtickformat('degrees')
ytickformat('degrees')
title('precipitation rate [m/yr] on day 1',...
    'FontWeight','Normal','FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
acc_plots
hold off
set(gca, 'nextplot','replacechildren', 'Visible','on');
vidObj = VideoWriter('movie_PRECIP_mean.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 32;
open(vidObj);
writeVideo(vidObj, getframe(gcf));
%% end initial plot for movie

%% time stepping
for ii=1:num_steps
    output_T_Time_Series(ii) = mean(mean(output_T));
    PRECIP_TEMP = fcn_precip(output_T,rho_air,rho_sea,Hq,output_Q,dt);
    output_PRECIP_Time_Series(ii) = mean(mean(PRECIP_TEMP));
    fprintf('mean precip: %g  mean q: %g  at step %g\n',output_PRECIP_Time_Series(ii),mean(mean(output_Q)),ii)
    if (mod(ii,7)==1)||(mod(ii,13)==1)
        [T_new,Q_new,PRECIP] = ode_fw_fe(output_T,SST,output_Q,Co,...
            rho_air,rho_sea,epsO,Ha,Hq,Lnu,So,...
            sigma,Crhoa,DY,DX,Coal,epsP,epsA,...
            nu,kappa,S,U,sea,qsatSST,dt,PRECIP);
    else
        [T_new,Q_new,PRECIP] = ode_fw_lf(output_T,T_old,SST,output_Q,Q_old,Co,...
            rho_air,rho_sea,epsO,Ha,Hq,Lnu,So,...
            sigma,Crhoa,DY,DX,Coal,epsP,epsA,...
            nu,kappa,S,U,sea,qsatSST,dt,PRECIP);
    end
    T_old = output_T;
    output_T = T_new;
    Q_old = output_Q;
    output_Q = Q_new;
    output_mean_PRECIP = ((ii-1)*output_mean_PRECIP + PRECIP)./ii;
    %% plot
    if (mod(ii,16)==1)&&(ii<num_steps/5)
        contourf(LON_GRID,LAT_GRID,output_mean_PRECIP,'LineStyle','none','LevelList',z3);
        colorbar('eastoutside');
        hold on
        contour(LON_GRID,LAT_GRID,0+land,'Color','k')
        caxis(coloraxis3)
        axis(coords)
        xtickformat('degrees')
        ytickformat('degrees')
        title(['precipitation rate [m/yr] on day ',num2str(1+floor(ii*dt/86400))],...
            'FontWeight','Normal','FontSize',16)
        xlabel('Longitude')
        ylabel('Latitude')
        acc_movie
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gcf));
    elseif (mod(ii,30)==1)
        contourf(LON_GRID,LAT_GRID,output_mean_PRECIP,'LineStyle','none','LevelList',z3);
        colorbar('eastoutside');
        hold on
        contour(LON_GRID,LAT_GRID,0+land,'Color','k')
        caxis(coloraxis3)
        axis(coords)
        xtickformat('degrees')
        ytickformat('degrees')
        title(['precipitation rate [m/yr] on day ',num2str(1+floor(ii*dt/86400))],...
            'FontWeight','Normal','FontSize',16)
        xlabel('Longitude')
        ylabel('Latitude')
        acc_movie
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gcf));
    end
    %% end plot
    
    %% check for instability
    if  sum(sum(isnan(output_T)))>1
        disp('scheme has blown up!')
        break
    elseif  sum(sum(isnan(output_Q)))>1
        disp('scheme has blown up!')
        break
        %     elseif min(min(Q))<0
        %         disp('negative humidity!')
        %         break
    end
    %% end check for instability
end
%% end time stepping
close(vidObj);

save output_mean_PRECIP output_mean_PRECIP
toc()