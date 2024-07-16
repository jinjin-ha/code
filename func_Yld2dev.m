%%
function [DFDS]=func_Yld2dev(M,A,S)

% Derivative
[L1p, L2p]=func_L(A);
[X1p,X2p]=func_lineartrans(L1p,L2p,S);
[XP1p,H1]=func_principal(X1p); [XP2p,H2]=func_principal(X2p);

P(3)=XP1p(1)-XP1p(2);
P(1)=2*XP2p(2)+XP2p(1);
P(2)=2*XP2p(1)+XP2p(2);

% FIRST DERIVATIVES WITH RESPECT TO PRINCIPAL STRESSES
DFL1(1)=M*( sign(P(3))*abs(P(3))^(M-1));    
DFL1(2)=M*(-sign(P(3))*abs(P(3))^(M-1));    
DFL2(1)=M*( sign(P(1))*abs(P(1))^(M-1)+2*sign(P(2))*abs(P(2))^(M-1));
DFL2(2)=M*(2*sign(P(1))*abs(P(1))^(M-1)+ sign(P(2))*abs(P(2))^(M-1));

% FIRST TRANSFORMATION
D1=sqrt(H1(3));
if D1>(10^(-10))
    DLX1(1,1)=(1+H1(2)/D1)/2;
    DLX1(1,2)=(1-H1(2)/D1)/2;
    DLX1(1,3)=2*X1p(6)/D1;
    DLX1(2,1)=(1-H1(2)/D1)/2;
    DLX1(2,2)=(1+H1(2)/D1)/2;
    DLX1(2,3)=-2*X1p(6)/D1;
    Y1(1)=DFL1(1)*DLX1(1,1)+DFL1(2)*DLX1(2,1);
    Y1(2)=DFL1(1)*DLX1(1,2)+DFL1(2)*DLX1(2,2);
    Y1(6)=DFL1(1)*DLX1(1,3)+DFL1(2)*DLX1(2,3);
else
    Y1(1)=DFL1(1);
    Y1(2)=DFL1(1);
    Y1(6)=0.0;
end

% SECOND TRANSFORMATION
D2=sqrt(H2(3));
if D2>(10^(-10))
    DLX2(1,1)=(1+H2(2)/D2)/2;
    DLX2(1,2)=(1-H2(2)/D2)/2;
    DLX2(1,3)=2*X2p(6)/D2;
    DLX2(2,1)=(1-H2(2)/D2)/2;
    DLX2(2,2)=(1+H2(2)/D2)/2;
    DLX2(2,3)=-2*X2p(6)/D2;
    Y2(1)=DFL2(1)*DLX2(1,1)+DFL2(2)*DLX2(2,1);
    Y2(2)=DFL2(1)*DLX2(1,2)+DFL2(2)*DLX2(2,2);
    Y2(6)=DFL2(1)*DLX2(1,3)+DFL2(2)*DLX2(2,3);
else
    Y2(1)=DFL2(1);
    Y2(2)=DFL2(1);
    Y2(6)=0.0;
end

DFDS=Y1*L1p+Y2*L2p;
DFDS(6)=DFDS(6)/2;  % Transform shear derivatives in strains
end
