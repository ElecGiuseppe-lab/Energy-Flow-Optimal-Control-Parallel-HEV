%% Longitudinal vehicle dynamic model and energy optimization combined cycle
% Si determinano i punti operativi coppia-velocitá angolare a partire dal
% ciclo di guida e, quindi, la potenza che dovrá essere fornita dal
% powertrain alle ruote, mediante il calcolo della dinamica di base
% del veicolo e l'implementazione della strategia di controllo ottimo.

close all

Parametri_Hyundai_Tucson;

% Carico cicli di guida
load FUDS
load NEDC
load SFUDS
load WLTC_full

vel_FUDS = FUDS(:,2);
vel_SFUDS = SFUDS(:,2);
vel_WLTC_full = WLTC_full(:,2);
vel_NEDC = NEDC(:,2);

% Ciclo di guida COMBINATO
vel_cycle = cat(1,vel_FUDS,vel_WLTC_full,vel_SFUDS,vel_FUDS,vel_FUDS,vel_WLTC_full);
t_cycle = (1:numel(vel_cycle))';
vel_cycle_comb = vel_cycle;
t_cycle_comb = t_cycle;

dt_cycle = 1; % time step
%%

% Accelerazione [m/s^2]
acc_cycle = zeros(size(vel_cycle));

for i=2:numel(vel_cycle)
    if vel_cycle(i) > vel_max
        acc_cycle(i) = 0;
    else
        acc_cycle(i) = (vel_cycle(i) - vel_cycle(i-1))/dt_cycle;
    end
end

% Distanza percorsa [m]
x_cycle = zeros(size(vel_cycle));

for i=1:numel(vel_cycle)-1
    x_cycle(i+1) = x_cycle(i) + dt_cycle*(vel_cycle(i) + vel_cycle(i+1))/2;
end

x_cycle_comb = x_cycle;

% Plot velocità, accelerazione e distanza percorsa
figure(1)
subplot(3,1,1)
plot(t_cycle,vel_cycle);
set(xlabel('Time/sec'),'Interpreter','latex');
set(ylabel('Velocit\`a [m/s]'),'Interpreter','latex');
subplot(3,1,2)
plot(t_cycle,acc_cycle);
set(xlabel('Time/sec'),'Interpreter','latex');
set(ylabel('a [m/$s^2$]'),'Interpreter','latex');
subplot(3,1,3)
plot(t_cycle,x_cycle);
set(xlabel('Time/sec'),'Interpreter','latex');
set(ylabel('Distanza percorsa [m]'),'Interpreter','latex');
set(gcf,'Position',[200 100 750 550]);
%%

% Forze resistive esterne agenti sul veicolo in movimento
Fr = mhu_r*m*g*cosd(alpha); % Forza di resistenza al rotolamento [N]
Fg = m*g*sind(alpha); % Forza gravitazionale [N] (Fg=0 essendo su
                      % superficie piana)
Fa = 0.5*A*Cd*rho*vel_cycle.^2; % Forza di resistenza aerodinamica [N]
Fi = m_eff*acc_cycle; % Forza inerziale [N]

% Forza di trazione totale [N]
Ftr = zeros(size(vel_cycle));

for i=1:numel(vel_cycle)
    if vel_cycle(i) == 0 && acc_cycle(i) == 0
        Ftr(i) = 0;
    else
        Ftr(i) = Fr + Fg + Fa(i) + Fi(i);
    end
end

% Coppia, velocitá angolare e potenza richiesta alle ruote
T_wheel = Ftr*r_wheel;
w_wheel = vel_cycle/r_wheel;
P_wheel = T_wheel.*w_wheel;

% Regola del cambio
[G,SelGear] = TrasmissionModel(i_fd,i_gb,n_gb,vel_cycle);

% Coppia, velocitá angolare e potenza powertrain (punti operativi)
T_pwt = zeros(size(vel_cycle));
w_pwt = zeros(size(vel_cycle));
P_pwt = zeros(size(vel_cycle));

for i=1:numel(vel_cycle)
    % T>0 [Nm] é la coppia che dovrá essere fornita dal powertrain alle
    % ruote in fase di trazione. Il segno di T_wheel come esponenziale
    % dell'efficienza meccanica consente di tenere conto sia dei casi di
    % trazione che di frenata
    T_pwt(i) = T_wheel(i)/(G(i)*(eta_mech)^sign(T_wheel(i)));
    w_pwt(i) = w_wheel(i)*G(i); % velocità angolare [rad/sec]
    P_pwt(i) = T_pwt(i)*w_pwt(i); % P>0 é la potenza che dovrá essere
                                  % fornita dal powertrain alle ruote [W]
end

w_pwt_rpm = w_pwt*60/(2*pi); % [rpm]

% Curva coppia-velocitá ICE
[w_limite_ICE,T_limite_ICE] = InternalCombustionEngine(t_cycle,Fr);
% Curve di efficenza ICE
[val,W_ICE,T_ICE,eta_ICE_map,eta_ICE] = EfficencyICE(w_limite_ICE,...
    T_max_ICE,w_pwt,T_pwt); 
% Curva coppia-velocitá EM e GEN
[w_limite_MOT,T_limite_MOT,P_limite_MOT] = ElectricMotor(t_cycle,Fr); 
% Curve di efficienza EM e GEN
[Val,W,T,eta_MOT_map,eta_GEN_map,eta_EM] = EfficencyElectricMotor...
    (T_max_EM,w_limite_MOT,w_pwt,T_pwt); 

w_MOT = w_pwt;
w_ICE = w_pwt;
P_ICE = zeros(size(t_cycle));
for i = 1:numel(t_cycle)
    if P_pwt(i) > 0
        P_ICE(i) = P_pwt(i)/eta_ICE(i);
    else
        P_ICE(i) = 0;
    end
end

P_EM = P_pwt./(eta_EM'.^sign(T_pwt));

% Plot curva limite coppia-velocitá ICE + curve di efficienza
figure(2)
plot(w_limite_ICE,T_limite_ICE,'k','LineWidth',2);
hold on
contour(W_ICE,T_ICE,eta_ICE_map,val,'LineWidth',1.2,'ShowText','on');
colormap(turbo);
colorbar();
set(xlabel('${\omega}_{ICE}$ [rpm]'),'Interpreter','latex');
set(ylabel('$T_{ICE}$ [Nm]'),'Interpreter','latex');
set(legend('ICE', 'Efficienza',Location='bestoutside'),'Interpreter',...
    'latex');

% Plot curva limite coppia-velocitá EM + curve di efficienza
figure(3)
plot(w_limite_MOT,T_limite_MOT,'k','LineWidth',2); % curva limite MOT
hold on
plot(w_limite_MOT,-T_limite_MOT,'k','LineWidth',2); % curva limite GEN
hold on
[C,H]=contour(W,T,eta_MOT_map,Val,'LineWidth',1.2);
clabel(C,H,'fontsize',8);
colormap(turbo);
hold on
[C1,H1]=contour(W,-T,eta_GEN_map,Val,'LineWidth',1.2);
clabel(C1,H1,'fontsize',8);
colormap(turbo);
colorbar();
set(xlabel('${\omega}_{EM}$ [rpm]'),'Interpreter','latex');
set(ylabel('$T_{EM}$ [Nm]'),'Interpreter','latex');
set(legend('MOT','GEN', 'Efficienza',Location='bestoutside'),...
    'Interpreter','latex');
%%

% Energia consumata in accelerazione [J]
[E_cycle_pos_ICE,fc_comb] = ForwardMotoring(t_cycle,P_ICE,Q_lhv);

% Energia recuperabile con la frenata rigenerativa in decelerazione [J]
E_cycle_neg_GEN = ForwardBraking(vel_cycle,P_EM);

% Plot energia consumata e rigenerata
figure(4)
yyaxis left
plot(t_cycle,E_cycle_pos_ICE/10^6,'k'); % espresso in [MJ]
hold on
plot(t_cycle,abs(E_cycle_neg_GEN)/10^6,'g');
set(ylabel('Energia [MJ]'),'Interpreter','latex');
yyaxis right
plot(t_cycle,fc_comb,'m');
set(xlabel('Time/sec'),'Interpreter','latex');
set(ylabel('Consumo di carburante [lt]'),'Interpreter','latex');
set(legend('Energia consumata solo ICE',...
    'Energia recuperata con frenata rigenerativa',...
    'Consumo di carburanre',Location='northwest'),'Interpreter','latex');
%%

% Modello circuitale di ordine 0th batteria Li-ion a temperatura ambiente (25 °C)
[SoC,Voc_comb,Voc_8th,R_int_d,R_int_c,Qp,I_batt] = BatteryModel(t_cycle);

% Plot tensione a vuoto e ristenza interna in funzione del SoC
figure(5)
yyaxis left
plot(SoC,Voc_8th);
set(ylabel('$V_{oc}$ [V]'),'Interpreter','latex');
yyaxis right
plot(SoC,R_int_c,'r');
hold on
plot(SoC,R_int_d,'r');
set(ylabel('$R_{int}$ [Ohm]'),'Interpreter','latex');
set(xlabel('SoC'),'Interpreter','latex');
set(legend('$V_{oc}$','$R_{int,chg}$','$R_{int,dis}$',Location='northwest'),'Interpreter','latex');

% Plot capacitá batteria in funzione della corrente di scarica
figure(6)
plot(I_batt,Qp);
set(ylabel('$Q_{batt}$ [Ah]'),'Interpreter','latex');
set(xlabel('$I_{batt}$ [A]'),'Interpreter','latex');
set(legend('Modified Peukert',Location='northeast'),'Interpreter','latex');
ylim([0 24]);
%%

% Controllo ottimo: Pontryagin's Minimum Principle (shooting method)
[SoC_end,vec_lambda,SoC_cycle] =PontryaginMinimumPrinciple(t_cycle,w_pwt,P_pwt,P_limite_MOT);

% Bisezione per ricavare il valore iniziale ottimale del co-stato per il quale é soddisfatta la condizione sul SoC target
[lambda_opt,SoC_opt,P_ICE_cycle_opt,var_u_opt,P_fuel_cycle_opt,eta_batt_cycle,lamb_cycle_opt] =Bisezione(SoC_end,vec_lambda,t_cycle,...
    w_pwt,P_pwt,P_limite_MOT); 
%%

% Equivalent consumption Minimization Strategy (ECMS)
[SoC_end_ECMS,SoC_ECMS,P_ICE_cycle_ECMS,var_u_ECMS,P_fuel_cycle_ECMS,lamb_cycle_ECMS]=ECMS(eta_batt_cycle,lambda_opt,t_cycle,...
    w_pwt,P_pwt,P_limite_MOT);

% Continuos Adaptive Equivalent consumption Minimization Strategy (CA-ECMS)
lambda_CAECMS = lambda_opt-0.2; % co-stato iniziale
Kp = 8;
Ki = 0.4;
[SoC_end_CAECMS,SoC_CAECMS,P_ICE_cycle_CAECMS,var_u_CAECMS,P_fuel_cycle_CAECMS,lamb_cycle_CAECMS] =ContinuousAdaptiveECMS(Kp,Ki,...
    lambda_CAECMS,t_cycle,w_pwt,P_pwt,P_limite_MOT);

% Discrete Adaptive Equivalent consumption Minimization Strategy (DA-ECMS)
lambda_DAECMS = lambda_opt-0.2; % co-stato iniziale (non ottimale)
iter_update = 500; % numero di iterazioni prima di aggiornare il co-stato
Kd = 2; % coefficiente di proporzionalitá
[SoC_end_DAECMS,SoC_DAECMS,P_ICE_cycle_DAECMS,var_u_DAECMS,P_fuel_cycle_DAECMS,lamb_cycle_DAECMS] =DiscreteAdaptiveECMS(iter_update,Kd,...
    lambda_DAECMS,t_cycle,w_pwt,P_pwt,P_limite_MOT);
%%

%Plot traiettorie SoC PMP
figure(7)
SoC_opt_comb = SoC_opt;
plot(t_cycle,SoC_cycle','b',t_cycle,SoC_opt_comb,'r');
hold on
yline(SoC_target,'--');
set(ylabel('SoC'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');

% Plot effetto del valore iniziale del co-stato sul SoC finale
figure(8)
errore_comb = SoC_cycle(:,end)-SoC_target;
plot(vec_lambda,errore_comb,'.-b',MarkerSize=11);
hold on
errore_opt_comb = SoC_opt(end)-SoC_target;
plot(lambda_opt,errore_opt_comb,'.r',MarkerSize=11);
hold on
xline(lambda_opt,'--r');
set(xlabel('${\lambda}_0$'),'Interpreter','latex');
set(ylabel('Errore (${\Delta}SoC$)'),'Interpreter','latex');
yline(0,'--');
text(vec_lambda(1,1)+0.01,0.02,'$SoC_{end}$ = $SoC_{target}$',...
    Interpreter='latex')
%%

% Plot traiettoria SoC ottimale PMP vs ECMS vs CA-ECMS vs DA-ECMS
figure(9)
plot(t_cycle,SoC_opt_comb,ColorMode="manual",Color="#A2142F");
hold on
SoC_ECMS_comb = SoC_ECMS;
plot(t_cycle,SoC_ECMS_comb,ColorMode="manual",Color="#77AC30");
hold on
SoC_CAECMS_comb = SoC_CAECMS;
plot(t_cycle,SoC_CAECMS_comb,ColorMode="manual",Color="#0072BD");
hold on
SoC_DAECMS_comb = SoC_DAECMS;
plot(t_cycle,SoC_DAECMS_comb,ColorMode="manual",Color="#D95319");
hold on
yline(SoC_target,'--');
text(100,SoC_target+0.002,'SoC target')
set(ylabel('SoC'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
set(legend('PMP','ECMS','Continuous A-ECMS','Discrete A-ECMS',...
    Location='bestoutside'),'Interpreter','latex');
set(gcf,'Position',[800 200 1200 500]);

%%

% Punti operativi EM e ICE (PMP, ECMS, CA-ECMS, DA-ECMS)
P_ICE_cycle_opt_comb = P_ICE_cycle_opt;
P_EM_cycle_opt_comb = var_u_opt;

P_ICE_cycle_comb_ECMS = P_ICE_cycle_ECMS;
P_EM_cycle_comb_ECMS = var_u_ECMS;

P_ICE_cycle_comb_CAECMS = P_ICE_cycle_CAECMS;
P_EM_cycle_comb_CAECMS = var_u_CAECMS;

P_ICE_cycle_comb_DAECMS = P_ICE_cycle_DAECMS;
P_EM_cycle_comb_DAECMS = var_u_DAECMS;

T_ICE_cycle_opt = P_ICE_cycle_opt./w_pwt'; % Punti operativi ICE (PMP)
T_ICE_cycle_opt(T_ICE_cycle_opt==0) = nan;
T_EM_cycle_opt = var_u_opt./w_pwt'; % Punti operativi MOT (PMP)
T_EM_cycle_opt(T_EM_cycle_opt==0) = nan;

% Plot punti operativi EM e ICE
figure(10)
plot(w_pwt_rpm,T_ICE_cycle_opt,'.');
hold on
plot(w_limite_ICE,T_limite_ICE,'k','LineWidth',2);
hold on
contour(W_ICE,T_ICE,eta_ICE_map,val,'LineWidth',1.2,'ShowText','on');
colormap(turbo);
colorbar();
set(xlabel('${\omega}_{ICE}$ [rpm]'),'Interpreter','latex');
set(ylabel('$T_{ICE}$ [Nm]'),'Interpreter','latex');

figure(11)
plot(w_pwt_rpm,T_EM_cycle_opt,'.');
hold on
plot(w_limite_MOT,T_limite_MOT,'k','LineWidth',2); % Curva limite EM
hold on
plot(w_limite_MOT,-T_limite_MOT,'k','LineWidth',2); % Curva limite GEN
hold on
[C,H]=contour(W,T,eta_MOT_map,Val,'LineWidth',1.2);
clabel(C,H,'fontsize',8);
colormap(turbo);
hold on
[C1,H1]=contour(W,-T,eta_GEN_map,Val,'LineWidth',1.2);
clabel(C1,H1,'fontsize',8);
colormap(turbo);
colorbar();
set(xlabel('${\omega}_{EM}$ [rpm]'),'Interpreter','latex');
set(ylabel('$T_{EM}$ [Nm]'),'Interpreter','latex');
%%

% Consumo di carburante PMP, ECMS, CA-ECMS
E_ice_cycle = zeros(1,numel(t_cycle));
E_ice_cycle_ECMS = zeros(1,numel(t_cycle));
E_ice_cycle_CAECMS = zeros(1,numel(t_cycle));
E_ice_cycle_DAECMS = zeros(1,numel(t_cycle));

for j=1:numel(t_cycle)-1
    E_ice_cycle(j+1) = E_ice_cycle(j) + dt_cycle*P_fuel_cycle_opt(j);
    E_ice_cycle_ECMS(j+1) = E_ice_cycle_ECMS(j) + dt_cycle*P_fuel_cycle_ECMS(j);
    E_ice_cycle_CAECMS(j+1) = E_ice_cycle_CAECMS(j) + dt_cycle*P_fuel_cycle_CAECMS(j);
    E_ice_cycle_DAECMS(j+1) = E_ice_cycle_DAECMS(j) + dt_cycle*P_fuel_cycle_DAECMS(j);
end

fc_comb_opt = E_ice_cycle/Q_lhv;
fc_comb_ECMS = E_ice_cycle_ECMS/Q_lhv;
fc_comb_CAECMS = E_ice_cycle_CAECMS/Q_lhv;
fc_comb_DAECMS = E_ice_cycle_DAECMS/Q_lhv;

% Risparmio di carburante al traguardo PMP, ECMS, CA-ECMS
fc_sav_comb = (fc_comb(end)-fc_comb_opt(end))/fc_comb(end)*100;
fc_sav_comb_ECMS = (fc_comb(end)-fc_comb_ECMS(end))/fc_comb(end)*100; 
fc_sav_comb_CAECMS = (fc_comb(end)-fc_comb_CAECMS(end))/fc_comb(end)*100; 
fc_sav_comb_DAECMS = (fc_comb(end)-fc_comb_DAECMS(end))/fc_comb(end)*100; 
