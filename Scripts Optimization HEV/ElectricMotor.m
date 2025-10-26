%% Modello motore elettrico
% Funzione che consente di ottenere la curva caratteristica coppia-velocitá del motore elettrico in funzione dei parametri di progetto

function [w_limite_MOT,T_limite_MOT,P_limite_MOT] =ElectricMotor(t_cycle,Fr)

Parametri_Hyundai_Tucson;

% Curva limite coppia motore [Nm]
vel1 = zeros(size(t_cycle));
dt = 1;
Ftr_limite = zeros(size(t_cycle));

for j=1:numel(t_cycle)-1
    % Prima fase di accelerazione a coppia costante (T = T_max)
        if vel1(j) <= vel_cr_EM
            % Forza di trazione massima a coppia costante [Nm]
            Ftr_limite(j) = T_max_EM*G_em/r_wheel*eta_mech; 
            vel1(j+1) = vel1(j) + dt*((Ftr_limite(j) - Fr)/m_eff-(0.5*A*Cd*rho/m_eff)*vel1(j)^2);
        
            % Seconda fase di accelerazione a potenza costante (P = P_max)
        elseif vel1(j) > vel_cr_EM && vel1(j) <= vel_max
            % Forza di trazione a potenza costante [Nm]
            Ftr_limite(j) = P_max_EM/vel1(j)*eta_mech; 
            vel1(j+1) = vel1(j) + dt*((Ftr_limite(j) - Fr)/m_eff-(0.5*A*Cd*rho/m_eff)*vel1(j)^2);
        elseif vel1(j) > vel_max
            vel1(j+1) = vel1(j);
        end
end

w_limite_MOT = vel1*(G_em/r_wheel)*60/(2*pi); % espressa in [rpm]

T_limite_MOT = zeros(size(t_cycle));
for i=1:numel(t_cycle)
    if vel1(i)<=vel_cr_EM
        T_limite_MOT(i)=T_max_EM;
    elseif vel1(i) > vel_cr_EM
        T_limite_MOT(i)=P_max_EM/(vel1(i)*G_em/r_wheel);
    end
end

P_limite_MOT = zeros(size(t_cycle));
for i=1:numel(t_cycle)
    if vel1(i)<=vel_cr_EM
        P_limite_MOT(i)=T_max_EM*G_em/r_wheel*vel1(i);
    elseif vel1(i) > vel_cr_EM
        P_limite_MOT(i)=T_limite_MOT(i)*G_em/r_wheel*vel1(i);
    end
end

end