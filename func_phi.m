%% SIGBAR
% YF: Yield function type: [n=1] von-Mises, [n=2] Hosford, [n=3] Yld2000-2d
% m: Yield function exponent
% a: Yield function parameter
% S=[S11 S22 S33 S23 S31 S12]
function [SB]=func_phi(YF,m,k,a,S)

if YF==1 % von-Mises 
    [XP,~]=func_3Dprin(S);
% SB=sqrt(((S(1)-S(2))^2+(S(2)-S(3))^2+(S(3)-S(1))^2+6*(S(4)^2+S(5)^2+S(6)^2)/2);
    SBm=((XP(1)-XP(2))^2+(XP(2)-XP(3))^2+(XP(3)-XP(1))^2)/2;
    SB=(SBm)^(1/2);
end

if YF==2 % Hosford
    [XP,~]=func_3Dprin(S);
    SBm=((XP(1)-XP(2))^m+(XP(2)-XP(3))^m+(XP(3)-XP(1))^m)/2;
    SB=(SBm)^(1/m);
end

if YF==3 % Yld2000-2d
    [L1p, L2p]=func_L(a);
    [X1p,X2p]=func_lineartrans(L1p,L2p,S);
%     [XP1p,~]=func_3Dprin(X1p); [XP2p,~]=func_3Dprin(X2p);
    [XP1p,~]=func_principal(X1p); [XP2p,~]=func_principal(X2p); % Plane-stress principal value is necessary!
    P(3)=XP1p(1)-XP1p(2);
    P(1)=2*XP2p(2)+XP2p(1);
    P(2)=2*XP2p(1)+XP2p(2);    
    SBm=((P(3))^m+(P(1))^m+(P(2))^m)/k;
    SB=(SBm)^(1/m);
end

if abs(SB)<1D-16; SB = 0; end
    
% Not implemented
% Yld2004-18p
% Hill48
end