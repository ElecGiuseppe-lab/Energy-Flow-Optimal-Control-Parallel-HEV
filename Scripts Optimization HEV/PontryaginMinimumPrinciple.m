%% Strategia di controllo ottimale: Pontryagin's minimum principle
% Minimizzazione della funzione Hamiltoniana (consumo carburante) con controllo sul SoC della batteria, soddisfando vincoli globali e locali.
% Un metodo standard per risolvere un problema di controllo ottimo utilizzando Pontryagin é il cosiddetto "Shooting Method"

function [SoC_end,vec_lambda,SoC_cycle] =PontryaginMinimumPrinciple(t_cycle,w_pwt,P_pwt,P_limite_MOT)

Parametri_Hyundai_Tucson;

dt =1; % time step

% Inizializzazione del co-stato
vec_lambda = 2.7:0.01:3.05;
lambda = (vec_lambda.*ones(size(t_cycle)))';

% Funzione di penalitá
w_factor = zeros(numel(vec_lambda),numel(t_cycle)); 

% Inizializzazione SoC
SoC_cycle = zeros(numel(vec_lambda),numel(t_cycle));
SoC_cycle(:,1) = SoC_initial;

% Inizializzazione energia immagazzinata nella batteria
E_ech_cycle = E_batt.*SoC_cycle;

% Inizializzazione tensione a vuoto
Voc_cycle = zeros(numel(vec_lambda),numel(t_cycle));
Voc_cycle(:,1) = 65*(k_8th(1)+k_8th(2)*SoC_cycle(:,1)+k_8th(3)*SoC_cycle(:,1).^2+k_8th(4)*SoC_cycle(:,1).^3+...
    +k_8th(5)*SoC_cycle(:,1).^4+k_8th(6)*SoC_cycle(:,1).^5+k_8th(7)*SoC_cycle(:,1).^6+k_8th(8)*SoC_cycle(:,1).^7+...
    +k_8th(9)*SoC_cycle(:,1).^8);

% Inizializzazione resistenza interna in fase di carica
Rint_ch_cycle = zeros(numel(vec_lambda),numel(t_cycle));
Rint_ch_cycle(:,1) = (b_ch(1)*SoC_cycle(:,1).^4+b_ch(2)*SoC_cycle(:,1).^3+b_ch(3)*SoC_cycle(:,1).^2+...
    +b_ch(4)*SoC_cycle(:,1)+b_ch(5))*b_ch(6)*exp(b_ch(7)/(Te-b_ch(8)));

% Inizializzazione resistenza interna in fase di scarica
Rint_di_cycle = zeros(numel(vec_lambda),numel(t_cycle));
Rint_di_cycle(:,1) = (b_di(1)*SoC_cycle(:,1).^4+b_di(2)*SoC_cycle(:,1).^3+b_di(3)*SoC_cycle(:,1).^2+...
    +b_di(4)*SoC_cycle(:,1)+b_di(5))*b_di(6)*exp(b_di(7)/(Te-b_di(8)));


% Ottimizzazione
for i = 1:numel(vec_lambda) % ciclo sui valori iniziali del co-stato
    for j = 2:numel(t_cycle) % ciclo sugli istanti del ciclo di guida 
        if SoC_cycle(i,j-1) < lb_SoC
            w_factor(i,j) = wSoC;    
        elseif SoC_cycle(i,j-1) > up_SoC 
            w_factor(i,j) = -wSoC;
        else
            w_factor(i,j)=0;
        end

        if P_pwt(j) == 0
            var_u_opt(i,j) = 0; % variabile di controllo ottimale (definisce la potenza erogata dall'EM)
            I_batt_cycle(i,j) = 0; % corrente istantanea batteria
            P_fuel_cycle(1,j) = 0; % potenza istantanea ICE reale (tenendo conto dell'efficienza)
            P_ICE_cycle(i,j) = 0; % potenza istantanea ICE
            P_ech_cycle(i,j) = 0; % potenza istantanea generata dalla batteria
            eta_ICE_cycle(i,j) = 0; % efficienza istantanea ICE
            eta_EM_cycle(i,j) = 0; % efficienza istantanea EM
            eta_batt_cycle(i,j) = 0; % efficienza istantanea batteria
            SoC_cycle(i,j) = SoC_cycle(i,j-1); % stato di carica istantaneo
            E_ech_cycle(i,j) = E_ech_cycle(i,j-1); % energia batteria
            Voc_cycle(i,j) = Voc_cycle(i,j-1); % tensione istantanea batteria
            Rint_di_cycle(i,j) = Rint_di_cycle(i,j-1); % resistenza interna istantanea batteria (scarica)
            Rint_ch_cycle(i,j) = Rint_ch_cycle(i,j-1); % resistenza interna istantanea batteria (carica)
            lambda(i,j) = lambda(i,j-1);

        elseif P_pwt(j) > 0
            var_u = linspace(0,min(P_limite_MOT(j),P_pwt(j)),40);
            T_pwt = var_u./w_pwt(j);

            eta_EM = (w_pwt(j).*abs(T_pwt))./(w_pwt(j).*abs(T_pwt)+kc*T_pwt.^2 + ki*w_pwt(j) + kw*w_pwt(j)^3 + c_EM);

            eta_ICE = 0.42-2*10^-6.*abs(w_pwt(j)-2500/60*2*pi).^1.95-0.02/50^2.2.*abs((P_pwt(j)-var_u)./w_pwt(j)-200).^2;

            P_fuel = (P_pwt(j)-var_u)./eta_ICE; % potenza ICE

            I_batt = Voc_cycle(i,j-1)/(2*Rint_di_cycle(i,j-1))-sqrt((Voc_cycle(i,j-1)/(2*Rint_di_cycle(i,j-1)))^2-...
                    +var_u/Rint_di_cycle(i,j-1));

            eta_batt = (Qm.*(1-I_batt./I1)./((1-I_batt./I1)+(I_batt./I0).^n))./Q_batt;

            P_ech = var_u./eta_EM./eta_batt;
            Ham = P_fuel+(lambda(i,j-1)+w_factor(i,j))*P_ech; % Calcolo Hamiltoniana

            % Minimo Hamiltoniana per individuare il valore ottimo della variabile di controllo
            [Ham_min(i,j), index_min] = min(Ham);

            % Valori ottimali delle variabili per cui l'Hamiltonia risulta minimizzata
            var_u_opt(i,j) = var_u(index_min);
            I_batt_cycle(i,j) = I_batt(index_min);
            P_fuel_cycle(i,j) = P_fuel(index_min);
            P_ICE_cycle(i,j) = P_pwt(j)-var_u_opt(i,j);
            P_ech_cycle(i,j) = P_ech(index_min);
            eta_ICE_cycle(i,j) = eta_ICE(index_min);
            eta_EM_cycle(i,j) = eta_EM(index_min);
            eta_batt_cycle(i,j) = eta_batt(index_min);
            E_ech_cycle(i,j) = E_ech_cycle(i,j-1)-P_ech_cycle(i,j)/3600*dt;
            SoC_cycle(i,j) = E_ech_cycle(i,j)/E_batt;

            Voc_cycle(i,j) = 65*(k_8th(1)+k_8th(2)*SoC_cycle(i,j)+k_8th(3)*SoC_cycle(i,j).^2+k_8th(4)*SoC_cycle(i,j).^3+...
                +k_8th(5)*SoC_cycle(i,j).^4+k_8th(6)*SoC_cycle(i,j).^5+k_8th(7)*SoC_cycle(i,j).^6+k_8th(8)*SoC_cycle(i,j).^7+...
                +k_8th(9)*SoC_cycle(i,j).^8);

            Rint_di_cycle(i,j) = (b_di(1)*SoC_cycle(i,j)^4+b_di(2)*SoC_cycle(i,j)^3+b_di(3)*SoC_cycle(i,j)^2+...
                +b_di(4)*SoC_cycle(i,j)+b_di(5))*b_di(6)*exp(b_di(7)/(Te-b_di(8)));

            Rint_ch_cycle(i,j) = (b_ch(1)*SoC_cycle(i,j)^4+...
                +b_ch(2)*SoC_cycle(i,j)^3+b_ch(3)*SoC_cycle(i,j)^2+...
                +b_ch(4)*SoC_cycle(i,j)+...
                +b_ch(5))*b_ch(6)*exp(b_ch(7)/(Te-b_ch(8)));

            dSoC_SoC = 5e-4/35e3*P_ech_cycle(i,j); % determina la variazione del co-stato
            lambda(i,j) = lambda(i,j-1)-dt*(lambda(i,j-1)+w_factor(i,j))*dSoC_SoC; % aggiornamento del co-stato
%             lambda(i,j) = lambda(i,j-1);

        else
            var_u = linspace(max(-P_limite_MOT(j),P_pwt(j)),0,40);
            T_pwt = var_u./w_pwt(j);
            eta_EM = max(0,(w_pwt(j).*abs(T_pwt)-kc*T_pwt.^2- ki*w_pwt(j)-kw*w_pwt(j)^3-c_GEN)./(w_pwt(j).*abs(T_pwt)));

            I_batt = Voc_cycle(i,j-1)/(2*Rint_ch_cycle(i,j-1))-sqrt((Voc_cycle(i,j-1)/(2*Rint_ch_cycle(i,j-1)))^2-...
                    +var_u/Rint_ch_cycle(i,j-1));

            eta_batt = (Qm.*(1-abs(I_batt)./I1)./((1-abs(I_batt)./I1)+(abs(I_batt)./I0).^n))./Q_batt;
            
            P_ech = var_u.*eta_EM.*eta_batt;
            P_fuel = 0;
            Ham = P_fuel+(lambda(i,j-1)+w_factor(i,j))*P_ech;

            [Ham_min(i,j), index_min] = min(Ham);

            var_u_opt(i,j) = var_u(index_min);
            I_batt_cycle(i,j) = I_batt(index_min);
            P_fuel_cycle(i,j) = 0;
            P_ICE_cycle(i,j) = 0;
            P_ech_cycle(i,j) = P_ech(index_min); 
            eta_ICE_cycle(i,j) = 0;
            eta_EM_cycle(i,j) = eta_EM(index_min);
            eta_batt_cycle(i,j) = eta_batt(index_min);
            E_ech_cycle(i,j) = E_ech_cycle(i,j-1)-P_ech_cycle(i,j)/3600*dt;
            SoC_cycle(i,j) = E_ech_cycle(i,j)/E_batt;

            Voc_cycle(i,j) = 65*(k_8th(1)+k_8th(2)*SoC_cycle(i,j)+k_8th(3)*SoC_cycle(i,j).^2+k_8th(4)*SoC_cycle(i,j).^3+...
                +k_8th(5)*SoC_cycle(i,j).^4+k_8th(6)*SoC_cycle(i,j).^5+k_8th(7)*SoC_cycle(i,j).^6+k_8th(8)*SoC_cycle(i,j).^7+...
                +k_8th(9)*SoC_cycle(i,j).^8);

            Rint_di_cycle(i,j) = (b_di(1)*SoC_cycle(i,j)^4+b_di(2)*SoC_cycle(i,j)^3+b_di(3)*SoC_cycle(i,j)^2+...
                +b_di(4)*SoC_cycle(i,j)+b_di(5))*b_di(6)*exp(b_di(7)/(Te-b_di(8)));

            Rint_ch_cycle(i,j) = (b_ch(1)*SoC_cycle(i,j)^4+b_ch(2)*SoC_cycle(i,j)^3+b_ch(3)*SoC_cycle(i,j)^2+...
                +b_ch(4)*SoC_cycle(i,j)+b_ch(5))*b_ch(6)*exp(b_ch(7)/(Te-b_ch(8)));

            dSoC_SoC = 5e-4/35e3*P_ech_cycle(i,j);
            lambda(i,j) = lambda(i,j-1)-dt*(lambda(i,j-1)+w_factor(i,j))*dSoC_SoC;
%             lambda(i,j) = lambda(i,j-1);
         end
    end
    % Valore finale di ogni traiettoria del SoC in funzione del co-stato iniziale
    SoC_end(i) = SoC_cycle(i,end);
end
end


