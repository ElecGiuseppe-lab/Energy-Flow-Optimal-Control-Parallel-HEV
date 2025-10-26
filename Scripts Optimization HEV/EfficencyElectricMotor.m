%% Mappa di efficienza MOT e GEN ed efficienza nei punti operativi

function [Val,W,T,eta_MOT_map,eta_GEN_map,eta_EM]=...
    EfficencyElectricMotor(T_max_EM,w_limite_MOT,w_pwt,T_pwt)

kc=0.1; % Coefficiente perdite rame
ki=0.5e-2; % Coefficiente perdite metallo
kw=1.2*10^-5; % Coefficiente perdite aerodinamiche
c_EM=20;
c_GEN=200;

W = linspace(min(w_limite_MOT),max(w_limite_MOT),100);
T = linspace(0,T_max_EM,100);
[W,T]= meshgrid(W,T); % Griglia di punti sul piano coppia-velocità angolare
W_rad = (W*2*pi)/60; % Velocità angolare espressa in [rad/sec]
Out_power = W_rad.*T;
Loss_power_EM= kc*T.^2 + ki*W_rad + kw*W_rad.^3 + c_EM;
Loss_power_GEN= kc*T.^2 + ki*W_rad + kw*W_rad.^3 + c_GEN;
eta_MOT_map = Out_power./(Out_power+Loss_power_EM);
eta_GEN_map = Out_power./(Out_power+Loss_power_GEN);
Val = [0.70, 0.85, 0.90, 0.92, 0.94, 0.95, 0.96];

% Calcolo efficienza MOT e GEN sui punti di lavoro
OutPower = w_pwt.*abs(T_pwt);
TotPowerEM = OutPower + kc*T_pwt.^2 + ki*w_pwt + kw*w_pwt.^3 + c_EM;
TotPowerGEN = OutPower - kc*T_pwt.^2 - ki*w_pwt - kw*w_pwt.^3 - c_GEN;

for i=1:numel(T_pwt)
    if T_pwt(i) >= 0
        % Efficienza motore elettrico
        eta_EM(i) = OutPower(i)/TotPowerEM(i); 
    else
        % Efficienza generatore
        eta_EM(i) = max([0,TotPowerGEN(i)/OutPower(i)]); 
    end
end

end


