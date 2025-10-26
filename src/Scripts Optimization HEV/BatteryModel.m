%% Modello batteria
% Viene utilizzato un modello circuitale equivalente 0th-order per tener
% conto delle perdite resistive (effetto Joule).
% Si determinano Voc e R_int in funzione del SoC attraverso l'utilizzo di
% formule empiriche ricavate sperimentalmente supposto T=25°C.

function [SoC,Voc_comb,Voc_8th,R_int_d,R_int_c,Qp,I_batt]=...
    BatteryModel(t_cycle)

Parametri_Hyundai_Tucson;

% Stato di carica batteria
SoC = linspace(0,1,length(t_cycle));

% Tensione a vuoto (8th order)
Voc_8th = 65*(k_8th(1)+k_8th(2)*SoC+k_8th(3)*SoC.^2+k_8th(4)*SoC.^3+...
    +k_8th(5)*SoC.^4+k_8th(6)*SoC.^5+k_8th(7)*SoC.^6+k_8th(8)*SoC.^7+...
    +k_8th(9)*SoC.^8);

% Tensione a vuoto (esponenziale + polinomio)
Voc_comb = 63*(a(1)*exp(a(2)*SoC)+a(3)*SoC+a(4)*SoC.^2+a(5)*SoC.^3+a(6)); 

% Resistenza interna in fase di carica e scarica
R_int_d = (b_di(1)*SoC.^4+b_di(2)*SoC.^3+b_di(3)*SoC.^2+b_di(4)*SoC+...
    +b_di(5))*b_di(6)*exp(b_di(7)/(Te-b_di(8)));

R_int_c = (b_ch(1)*SoC.^4+b_ch(2)*SoC.^3+b_ch(3)*SoC.^2+b_ch(4)*SoC+...
    b_ch(5))*b_ch(6)*exp(b_ch(7)/(Te-b_ch(8)));

% Modello Peukert modificato per tener conto degli effetti Coulombiani
I_batt = linspace(0,260,length(t_cycle));
for i=1:numel(I_batt)
    Qp(i) = Qm*(1-I_batt(i)/I1)/((1-I_batt(i)/I1)+(I_batt(i)/I0)^n);
end

