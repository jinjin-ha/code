%%
function [XP,H]=func_principal(X)
H(1) =  X(1) + X(2);
H(2) =  X(1) - X(2);
H(3) =  H(2)^2 + 4*X(6)^2;

XP(1) = (H(1) + sqrt(H(3)))/2;
XP(2) = (H(1) - sqrt(H(3)))/2;
XP(3) = -XP(1)-XP(2);

for i=1:3
    if abs(XP(i))<10^(-16)
        XP(i)=0.0;
    end
end
end