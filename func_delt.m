%% DELTA in [rad]
% E=[E11 E22 E33 E23 E31 E12]
function [delt]=func_delt(EE,YE)
[EEN,~]=func_norm(EE);
[YEN,~]=func_norm(YE);

[DB]=func_double(EEN,YEN);
delt=real(acos(DB));   % delta in [rad]
end