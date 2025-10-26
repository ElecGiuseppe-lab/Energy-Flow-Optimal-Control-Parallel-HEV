%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energia consumata solo ICE (Potenza > 0) in accelerazione [J]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E_cycle_pos_ICE,fc] = ForwardMotoring(t_cycle,P_ICE,Q_lhv)

E_cycle_pos_ICE = zeros(size(t_cycle));
dt=1;

% Energia consumata dal motore a combustione
for i=1:numel(t_cycle)-1
        E_cycle_pos_ICE(i+1) = E_cycle_pos_ICE(i) + dt*P_ICE(i);
end
fc = E_cycle_pos_ICE/Q_lhv; % Consumo di carburante [lt]

end