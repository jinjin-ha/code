%% COS_THETA
% S=[S11 S22 S33 S23 S31 S12];
function [cos]=func_cos(S)
SR=[1 0 0 0 0 0]; 

[SRD,~]=func_deviat(SR);
[SD,~]=func_deviat(S);

[SRDN,~]=func_norm(SRD);
[SDN,~]=func_norm(SD);

[cos]=func_double(SDN,SRDN);

if abs(cos)<1D-16; cos = 0; end
end
% end