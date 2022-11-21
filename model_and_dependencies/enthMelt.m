function [T,phi] = enthMelt(nonH,mixCut,H_to_T,phi_fun,TWater_fun)%(nonH,latHeat,rho_w,rho_i,c_pi,c_pw,DT,H_to_T)
    
% T = nan(Grid.p.Ny,Grid.p.Nx);
    % phi = nan(Grid.p.Ny,Grid.p.Nx);
    T = nan(length(nonH),1);
    phi = nan(length(nonH),1);
    
    %% create zones based on enthalpy
    iLog = nonH<=0;
    mLog = nonH>0 & nonH<mixCut;%((rho_w*latHeat)/(rho_i*c_pi*DT));
    wLog = nonH >= mixCut;%((rho_w*latHeat)/(rho_i*c_pi*DT));
    
    %% ice equations
%     T(iLog) = T_b + nonH(iLog)/(rho_i*c_pi);
    T(iLog) = H_to_T(nonH(iLog));%1+nonH(iLog);
    phi(iLog) = 0;
    
    %% mixture equations
%     T(mLog) = T_b;
%     phi(mLog) = nonH(mLog)/(rho_w*latHeat);
%     phi(mLog) = nonH(mLog) * (rho_i*c_pi*DT)/(rho_w*latHeat);
    T(mLog) = 1;
    phi(mLog) = phi_fun(nonH(mLog));
    
    %% water equations
%     T(wLog) = nonH(wLog)* (rho_i*c_pi)/(rho_w*c_pw) - latHeat/(DT*c_pw) + 1;
%     T(wLog) = T_b+(nonH(wLog)-rho_w*latHeat)/(rho_w*c_pw);
    T(wLog) = TWater_fun(nonH(wLog));
    phi(wLog) = 1;
end