function [dy] = hires_tempvol(t,y,PD)

dy(1) = PD.coolingrate(t);
dy(2) = PD.ASadditionrate(t);

dy = dy(:);