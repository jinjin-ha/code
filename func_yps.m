%% YIELD FUNCTION PLANE STRAIN
function [YPSG,YPSL]=func_yps(MM,KK,A,pang)

y0=[1.11 0.55];
if MM==2; y0=[1.1547 0.5775]; end
if MM==6; y0=[1.1167 0.5584]; end
if MM==8; y0=[1.0894 0.5447]; end

% Least-Square for Yld2k-2d
y=lsqnonlin(@(y) func_costPSY2k(y,A,MM,KK,pang),y0,[],[],[]);
% Out put stress tensor is in a local coordinate system
% Transformation in to a global coordinate system
T=[cosd(pang) sind(pang);-sind(pang) cosd(pang)]; PSTRy=T'*[y(1) 0; 0 y(2)]*T;

YPSL=[y(1) y(2)];
YPSG=[PSTRy(1,1) PSTRy(2,2) 0 0 0 PSTRy(1,2)];
end