%% Parametri vettura (Hyundai Tucson Full-Hybrid)


% Parametri fisici
m=1649+70*3; % massa veicolo con carico [kg]
m_eff=m+40; % massa efficace [kg] (tiene conto dell'accelerazione rotazionale)
r_wheel=0.35; % raggio ruota [m]
Cd=0.33; % coefficiente aerodinamico
mhu_r=0.015; % coefficiente resistenza rotolamento
A=2.6; % area veicolo[m^2]
rho=1.225; % densità aria [kg/m^3]
i_fd = 3.32; % rapporto di trasmissione final drive
n_gb = [1 2 3 4 5 6]; % Select gear (marce)
i_gb = [4.639 2.826 1.841 1.386 1 0.772]; % rapporti di trasmissione gearbox (1st,2nd,3rd,4th,5th,6th)
eta_mech=0.95; % efficienza di trasmissione meccanica (prodotto tra efficienza frizione/cambio e quella del differenziale)
g=9.81; % accelerazione gravità [m/s^2]
alpha=0; % angolo pendenza stradale [°]
%%

% Parametri ICE (motore a combustione interna)
G_ice = 3.5;
Q_lhv = 33e6; % fattore di conversione tra energia e carburante [J/l]
T_max_ICE= 265; % coppia massima ICE [Nm]
w_cr_ICE = 4500/60*2*pi; % pulsazione critica ICE [rad/s]
vel_cr_ICE = w_cr_ICE/G_ice*r_wheel;
P_max_ICE=T_max_ICE*w_cr_ICE; % potenza massima motore IC [W]
vel_max = 193/3.6; % velocità massima [m/s]
%%

% Parametri EM (motore elettrico)
kc=0.1; % coefficiente perdite rame
ki=0.5e-2; % coefficiente perdite metallo
kw=1.2*10^-5; % coefficiente perdite aerodinamiche
c_EM=20;
c_GEN=200;
G_em = 7;
T_max_EM = 264; % coppia massima motore elettrico [Nm]
P_max_EM = 44.2e3; % potenza massima motore EM [W]
w_cr_EM = P_max_EM/T_max_EM; % pulsazione critica EM [rad/s] sotto la quale la coppia è costante e sopra la quale la potenza è costante
vel_cr_EM = w_cr_EM/G_em*r_wheel; % velocità critica [m/s]
%%

% Parametri batteria Li-ion e modello Peukert modificato
I0= 200;
I1 = 260;
n = 3.3;
a = [-0.402, -50.58, 0.8849, -1.662, 1.482, 3.574]; % coefficienti tensione a vuoto (esponenziale + polinomio)
k_8th = [2.85, 4.80, -17.8, 38.59, 4.91, -210.52, 433.98,-368.45, 115.81]; % coefficienti tensione a vuoto (8th order)
b_di = [1.298e-1, -2.892e-1, 2.273e-1, -7.216e-2, 8.980e-2, 7.613e-1,10.14, 2.608e2]; % coefficienti R_int (fase di scarica)
b_ch = [1.369e-1, -2.518e-1, 1.609e-1, -4.100e-2, 8.210e-2, 7.192e-1,33.91, 1.999e2]; % coefficienti R_int (fase di carica)
Te = 298.15; %temperatura [K]
eta_coul = 0.95; % efficienza batteria
V_oc = 270; % tensione a circuito aperto [V]
E_batt = 1.49e3*4; % energia immagazzinata [Wh]
Q_batt = E_batt/V_oc; % Capacitá nominale [Ah]
Qm = Q_batt*0.98; % capacitá effettiva [Ah]
%%

% Vincoli SoC
lb_SoC = 0.4; % limite inferiore SoC
up_SoC = 0.8; % limite superiore SoC
SoC_initial = 0.7; % stato carica iniziale
SoC_target = 0.7; % stato di carica finale desiderato
wSoC = 10^2.03; % funzione di penalitá