%% FULL(3D) PRINCIPAL COMPONENT
% S=[S11 S22 S33 S23 S31 S12];
function [RKS,H]=func_3Dprin(X)
% H1, H2, H3 CALCULATION
H(1) = X(1)+X(2)+X(3);
H(2) = X(1)*X(2)+X(2)*X(3)+X(3)*X(1)-X(4)^2-X(5)^2-X(6)^2;
H(3) = X(1)*X(2)*X(3)+2*X(4)*X(5)*X(6) - X(1)*X(4)^2-X(2)*X(5)^2-X(3)*X(6)^2;

% CALCULATION OF PRINCIPAL VALUES
p = [1 -H(1) H(2) -H(3)];
RKS = roots(p);

for I=1:3
    if abs(RKS(I))<1D-16
        RKS(I) = 0;
    end
end

end