%% COST FUNCTION FOR PLANE STRAIN
function FP=func_costPSY2k(y,A,MM,KK,pang)

% LINEAR TRANSFORMATION TENSOR
[Lp, Ld]=func_L(A);
LpPS=[Lp(1,1) Lp(1,2) 0; Lp(2,1) Lp(2,2) 0; 0 0 Lp(6,6)];
LdPS=[Ld(1,1) Ld(1,2) 0; Ld(2,1) Ld(2,2) 0; 0 0 Ld(6,6)];

% Stress
Rot=[cosd(pang) sind(pang); -sind(pang) cosd(pang)];    % Local to Global
RotP=Rot.';  % Transpose
SigG=RotP*[1 0; 0 y(2)/y(1)]*Rot;   % Stress in global: always uniaxial in local

Xp=LpPS*[SigG(1,1) SigG(2,2) SigG(1,2)]';
Xd=LdPS*[SigG(1,1) SigG(2,2) SigG(1,2)]';

Xp(6)=Xp(3);Xd(6)=Xd(3);    % Extend vector
[XP,~]=func_principal(Xp);
[XD,~]=func_principal(Xd);

phi(1)=abs(XP(1)-XP(2));
phi(2)=abs(2*XD(2)+XD(1));
phi(3)=abs(2*XD(1)+XD(2));

% Cost function for stress
FP(1)=(phi(1)^MM+phi(2)^MM+phi(3)^MM)^(1/MM)-KK^(1/MM)/y(1);

% r-value
SigGv=[SigG(1,1) SigG(2,2) 0 0 0 SigG(1,2)];    % Extend vector
[DFDS]=func_Yld2dev(MM,A,SigGv);
EpsG=[DFDS(1) DFDS(6); DFDS(6) DFDS(2)];
EpsL=Rot*EpsG*RotP;

% Cost function for r-value
FP(2)=EpsL(2,2);

end