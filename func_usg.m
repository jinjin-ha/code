%% STRESS TENSOR
% S=[S11 S22 S33 S23 S31 S12]
function [S]=func_usg(Ang,sr)   % usg=uniaxial stress in global coordinate
% 2X2 Stress tentor in local coordinate
SL=[sr 0; 0 0];

% Rotational tensor; Ang in [degree]
Rot=[cosd(Ang) sind(Ang); -sind(Ang) cosd(Ang)];    % Local to Global
RotP=Rot.';  % Transpose

% 6 Stress components in glocal coordinate
SG=RotP*SL*Rot;   % Stress in global: always uniaxial in local
S=[SG(1,1) SG(2,2) 0 0 0 SG(1,2)];
end