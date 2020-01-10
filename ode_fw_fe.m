function [Tnew,Qnew,PRECIP] = ode_fw_fe(T,SST,Q,Co,...
    rho_air,rho_sea,epsO,Ha,Hq,Lnu,So,...
    sigma,Crhoa,DY,DX,Coal,epsP,epsA,...
    nu,kappa,S,U,sea,qsatSST,dt,PRECIP)
%% forward Euler time integrator for FW model

%% forcing terms
Ch = fcn_Ch(T,SST,U);
QSSW = fcn_QSSW(So,S,Coal,Co);
QSH = sea.*fcn_QSH(SST,T,rho_air,Ch,Crhoa,U);
QRR = sea.*fcn_QRR(SST,T,epsO,epsA,sigma);
QLW = fcn_QLW(T,epsP,sigma);
% QLH = sea.*fcn_QLH(T,Lnu,rho_sea,rho_air,Hq,Q,dt);
QLH = fcn_QLH(T,Lnu,rho_sea,rho_air,Hq,Q,dt);
EVAP = fcn_evap(T,Q,rho_air,rho_sea,SST,U,qsatSST);
%% end forcing terms

%% FE step for forcing terms
coeff1 = 1/rho_air/Ha/Crhoa;
Tnew = T + coeff1*dt.*(QSSW - QLW + QRR +...
   QSH + QLH);
coeff2 = 1/rho_air/Hq;
Qnew = Q + coeff2*dt.*rho_sea.*(sea.*EVAP-PRECIP)./31536000;
%% end FE step

%% in between pole averaging
np = 2;
sp = 45;
Tnew(1,:) = mean(Tnew(np,:));
Tnew(46,:) = mean(Tnew(sp,:));
Qnew(1,:) = mean(Qnew(np,:));
Qnew(46,:) = mean(Qnew(sp,:));
%% end in between pole averaging

%% Matsuno for QT
order = 2;
QT = fcn_QT(Tnew,rho_air,Ha,Crhoa,nu,DX,DY,order);
Tstar = Tnew + dt*coeff1.*QT;
Tstar(1,:) = mean(Tstar(np,:));
Tstar(46,:) = mean(Tstar(sp,:));
QT = fcn_QT(Tstar,rho_air,Ha,Crhoa,nu,DX,DY,order);
Tnew = Tnew + dt*coeff1.*QT;
%% end Matsuno for QT

%% Matsuno for MT & Q
MT = fcn_MT(Qnew,rho_air,Hq,kappa,DX,DY,order);
Qstar = Qnew + coeff2*dt.*MT;
Qstar(1,:) = mean(Qstar(np,:));
Qstar(46,:) = mean(Qstar(sp,:));
MT = fcn_MT(Qstar,rho_air,Hq,kappa,DX,DY,order);
Qnew = Qnew + coeff2*dt.*MT;
%% end Matsuno for MT & Q

%% final pole averaging
Tnew(1,:) = mean(Tnew(np,:));
Tnew(46,:) = mean(Tnew(sp,:));
Qnew(1,:) = mean(Qnew(np,:));
Qnew(46,:) = mean(Qnew(sp,:));
%% end final pole averaging

%% fix humidity to reflect any precipitation
[PRECIP,Qnew] = fcn_precip(Tnew,rho_air,rho_sea,Hq,Qnew,dt);
%% end fix humidity

end

