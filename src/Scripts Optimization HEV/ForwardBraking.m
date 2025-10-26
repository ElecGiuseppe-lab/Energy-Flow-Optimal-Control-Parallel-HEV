%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energia recuperabile con la frenata rigenerativa (Potenza < 0) in decelerazione [J]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E_cycle_neg_GEN = ForwardBraking(t_cycle,P_EM)

E_cycle_neg_GEN = zeros(size(t_cycle));
dt=1;
for i=1:numel(t_cycle)-1
    if P_EM(i) < 0
        E_cycle_neg_GEN(i+1) = E_cycle_neg_GEN(i) + dt*P_EM(i);
    else
        E_cycle_neg_GEN(i+1) = E_cycle_neg_GEN(i);
    end
end
end