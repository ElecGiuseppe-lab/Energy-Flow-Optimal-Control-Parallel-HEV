%% Mappa di efficienza ICE ed efficienza nei punti operativi

function [val,W_ICE,T_ICE,eta_ICE_map,eta_ICE]=...
    EfficencyICE(w_limite_ICE,T_max_ICE,w_pwt,T_pwt)

W_ICE = linspace(min(w_limite_ICE),max(w_limite_ICE),100);
T_ICE = linspace(0,T_max_ICE,100);
[W_ICE,T_ICE]= meshgrid(W_ICE,T_ICE); % Griglia di punti sul piano
                                      % coppia-velocità angolare

W_rad = (W_ICE*2*pi)/60; % Velocità angolare espressa in [rad/sec]

eta_ICE_map = 0.42-2*10^-6.*abs(W_rad-2500/60*2*pi).^1.95-...
    +0.02/50^2.2.*abs(T_ICE-200).^2; 

val = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.42,0.44];

% Calcolo efficienza ICE sui punti di lavoro
eta_ICE = 0.42-2*10^-6.*abs(w_pwt-2300/60*2*pi).^1.95-...
    +0.02/50^2.2.*abs(T_pwt-200).^2;

end



