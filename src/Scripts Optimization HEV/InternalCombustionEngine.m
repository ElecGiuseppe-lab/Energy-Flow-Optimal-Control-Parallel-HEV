%% Modello motore a combustione interna
% Funzione che consente di ottenere la curva caratteristica coppia-velocitá dell'ICE in funzione dei parametri di progetto

function [w_limite_ICE,T_limite_ICE]=InternalCombustionEngine(t_cycle,Fr)

Parametri_Hyundai_Tucson;

% Curva limite coppia motore [Nm] 
vel1 = zeros(size(t_cycle));
dt = 1;

for j=1:numel(t_cycle)-1

    %Prima fase di accelerazione a coppia costante
        if vel1(j) <= vel_cr_ICE
            % Forza di trazione massima a coppia costante [Nm]
            Ftr_limite = T_max_ICE*G_ice/r_wheel*eta_mech; 
            vel1(j+1) = vel1(j) + dt*((Ftr_limite - Fr)/m_eff-...
                +(0.5*A*Cd*rho/m_eff)*vel1(j)^2);

        %Seconda fase di accelerazione a potenza costante
        elseif vel1(j) > vel_cr_EM && vel1(j) <= vel_max
            % Forza di trazione a potenza costante [Nm]
            Ftr_limite = P_max_ICE/vel1(j)*eta_mech; 
            vel1(j+1) = vel1(j) + dt*(Ftr_limite/m_eff - Fr/m_eff-...
                +(0.5*A*Cd*rho/m_eff)*vel1(j)^2);
        elseif vel1(j) > vel_max
            vel1(j+1) = vel1(j);
        end
end

w_limite_ICE = vel1*(G_ice/r_wheel)*60/(2*pi); % espressa in [rpm]

T_limite_ICE = zeros(size(t_cycle));
for i=1:numel(t_cycle)
    if vel1(i) <= vel_cr_ICE
        T_limite_ICE(i)=T_max_ICE;
    elseif vel1(i) > vel_cr_ICE
        T_limite_ICE(i)=P_max_ICE/(vel1(i)*G_ice/r_wheel);
    end
end

end