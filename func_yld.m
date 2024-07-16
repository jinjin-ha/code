%% YIELD FUNCTION R-VALUE & PLANE-STRAIN STRESS
function [Ysr,Yr,YPSG,YPSL]=func_yld(YF,M,K,A,Ang,PAng)
[~,j]=size(Ang);[~,i]=size(PAng);
if YF==1 % Isotropic: von-Mises
    Yr=[1 1 1 1 1 1 1 1]; Ysr=[1 1 1 1 1 1 1 1];
    YPSG=[1.1547 0.5774 0;          % PS00
          0.8661 0.8661 0.2886;     % PS45
          0.5774 1.1547 0];         % PS90
    YPSL=[1.1547 0.5774;1.1547 0.5774;1.1547 0.5774];
end

if YF==2 % Isotropic: Hosford
    Yr=[1 1 1 1 1 1 1 1]; Ysr=[1 1 1 1 1 1 1 1];
    if M==6
    YPSG=[1.1167 0.5584 0 0 0 0;          % PS00
          0.8376 0.8376 0 0 0 0.2792;     % PS45
          0.5584 1.1167 0 0 0 0];         % PS90
    YPSL=[1.1167 0.5584;1.1167 0.5584;1.1167 0.5584];
    end
    if M==8
    YPSG=[1.0894 0.5447 0 0 0 0;          % PS00
          0.8171 0.8171 0 0 0 0.2724;     % PS45
          0.5447 1.0894 0 0 0 0];         % PS90
    YPSL=[1.0894 0.5447;1.0894 0.5447;1.0894 0.5447];
    end
end

if YF==3    % Yld2k-2d
    
    for n=1:j       % Unixail tension
        % Tensor transformation: change code in the future
        S(1)=cosd(Ang(n))^2; 
        S(2)=sind(Ang(n))^2;
        S(3)=0.0;
        S(4)=0.0;
        S(5)=0.0;
        S(6)=cosd(Ang(n))*sind(Ang(n));
        [L1p, L2p]=func_L(A);
        [X1p,X2p]=func_lineartrans(L1p,L2p,S);
        [XP1p,~]=func_principal(X1p); [XP2p,~]=func_principal(X2p);
        F=abs(XP1p(1)-XP1p(2))^M + abs(2*XP2p(1)+XP2p(2))^M + abs(2*XP2p(2)+XP2p(1))^M;
        CT=(K/F)^(1/M);
        S=S*CT;
        [DFDS]=func_Yld2dev(M,A,S);
    
        % Tensor transformation: change code in the future
        Ysr(n)=S(1)*(cosd(Ang(n)))^2+S(2)*(sind(Ang(n)))^2+2*S(6)*sind(Ang(n))*cosd(Ang(n));
        Yr(n)=(DFDS(1)*(sind(Ang(n)))^2+DFDS(2)*(cosd(Ang(n)))^2-2*DFDS(6)*sind(Ang(n))*cosd(Ang(n)))/(-DFDS(1)-DFDS(2));        
    end
    
        n=j+1;       % Biaxial tension
        % Tensor transformation: change code in the future
        S=[1 1 0 0 0 0];
        [X1p,X2p]=func_lineartrans(L1p,L2p,S);
        [XP1p,~]=func_principal(X1p); [XP2p,~]=func_principal(X2p);
        F=abs(XP1p(1)-XP1p(2))^M + abs(2*XP2p(1)+XP2p(2))^M + abs(2*XP2p(2)+XP2p(1))^M;
        CT=(K/F)^(1/M);
        S=S*CT;
        [DFDS]=func_Yld2dev(M,A,S);
    
        % Tensor transformation: change code in the future
        Ysr(n)=CT;
        Yr(n)=DFDS(2)/DFDS(1);
        
    for n=1:i       % Plane strain
        [YPSGn,YPSLn]=func_yps(M,K,A,PAng(n));
        YPSG(n,:)=YPSGn;
        YPSL(n,:)=YPSLn;
        
    end
    
end
end
