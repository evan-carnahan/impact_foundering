function [visc] = barrViscPhase(nonT,phi,viscFun)

visc = nan(size(nonT));
if all(phi == 0)
    phi = zeros(size(nonT));
end
wLog = nonT > 1;
iLog = ~wLog;
visc(iLog) = viscFun(nonT(iLog),phi(iLog));
visc(wLog) = viscFun(1,1);