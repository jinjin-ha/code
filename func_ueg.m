%% STRAIN TENSOR
% E=[E11 E22 E33 E23 E31 E12]
function [EG]=func_ueg(Ang,r)   % ueg=uniaxial strain in global coordinate
% 2X2 Strain tentor in local coordinate
DE1=1;
DE2=-r/(1+r);
EE=[DE1 0; 0 DE2];

% Rotational tensor; Ang in [degree]
Rot=[cosd(Ang) sind(Ang); -sind(Ang) cosd(Ang)];    % Local to Global
RotP=Rot.';  % Transpose

% 6 Strain components in glocal coordinate
DE=RotP*EE*Rot;   % Stress in global: always uniaxial in local
EG=[DE(1,1) DE(2,2) -DE(1,1)-DE(2,2) 0 0 DE(1,2)];
end