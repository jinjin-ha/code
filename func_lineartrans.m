%% LINEAR TRANSFORMATION
% L: 6x6 Tensor
% S=[S11 S22 S33 S23 S31 S12];
function [X1p,X2p]=func_lineartrans(L1p,L2p,S)
X1p=L1p*S';
X2p=L2p*S';
end