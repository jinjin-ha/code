%% NORMALIZATION
% S=[S11 S22 S33 S23 S31 S12];
function [SN,N]=func_norm(S)
N=sqrt((S(1)^2+S(2)^2+S(3)^2+2*(S(4)^2+S(5)^2+S(6)^2)));
SN=S/N;
end