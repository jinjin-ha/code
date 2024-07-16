%% DOUBLE DOT
% S=[S11 S22 S33 S23 S31 S12];
function [C]=func_double(A,B)
% Bt=B.';         % Transpose
% C=trace(A*Bt);
C=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)+2*(A(4)*B(4)+A(5)*B(5)+A(6)*B(6));
end