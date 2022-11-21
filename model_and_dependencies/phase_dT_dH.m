function [dT_dH] = phase_dT_dH(nonH,mixCut,dT_dH_ice,volCpRat)%(nonH,latHeat,rho_w,rho_i,c_pi,c_pw,DT,c_p)

dT_dH = nan(length(nonH),1);
%% create zones based on enthalpy
iLog = nonH<=0;
mLog = nonH>0 & nonH<mixCut;%((rho_w*latHeat)/(rho_i*c_pi*DT));
wLog = nonH >= mixCut;%((rho_w*latHeat)/(rho_i*c_pi*DT));

% all ice location
% if strcmp(c_p,'cons')
%     dT_dH(iLog) = 1;
% elseif strcmp(c_p,'var')
dT_dH(iLog) = dT_dH_ice(nonH(iLog));%1./sqrt(1+(2*b*DT*nonH(iLog))./c_pb);
% end

% mixture region
dT_dH(mLog) = 0;
% as upper bound check 
% dT_dH(mLog) = dT_dH_ice(nonH(mLog));

% all water region
dT_dH(wLog) = volCpRat;%(rho_i*c_pi)/(rho_w*c_pw);









