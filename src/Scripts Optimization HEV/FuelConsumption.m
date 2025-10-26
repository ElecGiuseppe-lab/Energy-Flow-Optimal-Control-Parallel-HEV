%% Fuel consumption
% Consumo di carburante (pre-ottimizzazione e post-ottimizzazione) per differenti profili di guida (urbano, extra-urbano, combinato)

clc
clear all
close all

Main_ciclo_urbano;
Main_ciclo_Extraurbano;
Main_ciclo_combinato;
%%

Y_axis_1 = [round(fc_urb(end),3) round(fc_urb_opt(end),3) round(fc_urb_ECMS(end),3) round(fc_urb_CAECMS(end),3) round(fc_urb_DAECMS(end),3);
    round(fc_extraUrb(end),3) round(fc_extraUrb_opt(end),3) round(fc_extraUrb_ECMS(end),3) round(fc_extraUrb_CAECMS(end),3) round(fc_extraUrb_DAECMS(end),3);
    round(fc_comb(end),3) round(fc_comb_opt(end),3) round(fc_comb_ECMS(end),3) round(fc_comb_CAECMS(end),3) round(fc_comb_DAECMS(end),3)];

Y_axis_2 = [round(fc_sav_urb,3) round(fc_sav_urb_ECMS,3) round(fc_sav_urb_CAECMS,3) round(fc_sav_urb_DAECMS,3);
    round(fc_sav_extraUrb,3) round(fc_sav_extraUrb_ECMS,3) round(fc_sav_extraUrb_CAECMS,3) round(fc_sav_extraUrb_DAECMS,3);
    round(fc_sav_comb,3) round(fc_sav_comb_ECMS,3) round(fc_sav_comb_CAECMS,3) round(fc_sav_comb_DAECMS,3)];

X_axis = [1 2 3];

% SoC, potenza ICE e potenza EM (PMP)
figure
subplot(3,1,1)
yyaxis left
plot(t_cycle_urb,P_ICE_cycle_opt_urb/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_urb,P_EM_cycle_opt_urb/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_urb,SoC_opt_urb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('PMP ciclo urbano','FontWeight','bold');

subplot(3,1,2)
yyaxis left
plot(t_cycle_extraUrb,P_ICE_cycle_opt_extraUrb/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_extraUrb,P_EM_cycle_opt_extraUrb/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_extraUrb,SoC_opt_extraUrb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('PMP ciclo extra-urbano','FontWeight','bold');

subplot(3,1,3)
yyaxis left
plot(t_cycle_comb,P_ICE_cycle_opt_comb/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_comb,P_EM_cycle_opt_comb/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_comb,SoC_opt_comb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('PMP ciclo combinato','FontWeight','bold');
set(gcf,'Position',[200 100 1200 700]);

% SoC, potenza ICE e potenza EM (ECMS)
figure
subplot(3,1,1)
yyaxis left
plot(t_cycle_urb,P_ICE_cycle_urb_ECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_urb,P_EM_cycle_urb_ECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_urb,SoC_ECMS_urb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('ECMS ciclo urbano','FontWeight','bold');

subplot(3,1,2)
yyaxis left
plot(t_cycle_extraUrb,P_ICE_cycle_extraUrb_ECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_extraUrb,P_EM_cycle_extraUrb_ECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_extraUrb,SoC_ECMS_extraUrb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('ECMS ciclo extra-urbano','FontWeight','bold');

subplot(3,1,3)
yyaxis left
plot(t_cycle_comb,P_ICE_cycle_comb_ECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_comb,P_EM_cycle_comb_ECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_comb,SoC_ECMS_comb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('ECMS ciclo combinato','FontWeight','bold');
set(gcf,'Position',[200 100 1200 700]);

% SoC, potenza ICE e potenza EM (Continuous A-ECMS)
figure
subplot(3,1,1)
yyaxis left
plot(t_cycle_urb,P_ICE_cycle_urb_CAECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_urb,P_EM_cycle_urb_CAECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_urb,SoC_CAECMS_urb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('Continuous A-ECMS ciclo urbano','FontWeight','bold');

subplot(3,1,2)
yyaxis left
plot(t_cycle_extraUrb,P_ICE_cycle_extraUrb_CAECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_extraUrb,P_EM_cycle_extraUrb_CAECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_extraUrb,SoC_CAECMS_extraUrb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('Continuous A-ECMS ciclo extra-urbano','FontWeight','bold');

subplot(3,1,3)
yyaxis left
plot(t_cycle_comb,P_ICE_cycle_comb_CAECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_comb,P_EM_cycle_comb_CAECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_comb,SoC_CAECMS_comb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('Continuous A-ECMS ciclo combinato','FontWeight','bold');
set(gcf,'Position',[200 100 1200 700]);

% SoC, potenza ICE e potenza EM (Discrete A-ECMS)
figure
subplot(3,1,1)
yyaxis left
plot(t_cycle_urb,P_ICE_cycle_urb_DAECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_urb,P_EM_cycle_urb_DAECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_urb,SoC_DAECMS_urb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('Discrete A-ECMS ciclo urbano','FontWeight','bold');

subplot(3,1,2)
yyaxis left
plot(t_cycle_extraUrb,P_ICE_cycle_extraUrb_DAECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_extraUrb,P_EM_cycle_extraUrb_DAECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_extraUrb,SoC_DAECMS_extraUrb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('Discrete A-ECMS ciclo extra-urbano','FontWeight','bold');

subplot(3,1,3)
yyaxis left
plot(t_cycle_comb,P_ICE_cycle_comb_DAECMS/10^3,'Color', [0.9290 0.6940 0.1250 0.3]);
hold on
plot(t_cycle_comb,P_EM_cycle_comb_DAECMS/10^3,'Color', [0.4940 0.1840 0.5560 0.3]);
set(ylabel('Power [kW]'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
yyaxis right
plot(t_cycle_comb,SoC_DAECMS_comb,LineWidth=1.5,ColorMode="manual",Color="#77AC30");
set(ylabel('SoC'),'Interpreter','latex');
set(legend('ICE','EM','SoC',Location='bestoutside'),'Interpreter','latex');
title('Discrete A-ECMS ciclo combinato','FontWeight','bold');
set(gcf,'Position',[200 100 1200 700]);
%%

% SoC PMP, ECMS, CA-ECMS, DA-ECMS
figure
subplot(3,1,1)
yyaxis right
plot(t_cycle_urb,vel_cycle_urb,'Color', [0 0 0 0.3]);
set(ylabel('Velocit\`a [m/s]'),'Interpreter','latex');
yyaxis left
plot(t_cycle_urb,SoC_opt_urb,ColorMode="manual",Color="#A2142F");
hold on
plot(t_cycle_urb,SoC_ECMS_urb,ColorMode="manual",Color="#77AC30",LineStyle="-");
hold on
plot(t_cycle_urb,SoC_CAECMS_urb,ColorMode="manual",Color="#0072BD",LineStyle="-");
hold on
plot(t_cycle_urb,SoC_DAECMS_urb,ColorMode="manual",Color="#D95319",LineStyle="-");
hold on
yline(SoC_target,'--');
text(100,SoC_target+0.002,'SoC target')
set(ylabel('SoC'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
set(legend('PMP','ECMS','Continuous A-ECMS','Discrete A-ECMS',Location='bestoutside'),'Interpreter','latex');
title('Ciclo urbano','FontWeight','bold');

subplot(3,1,2)
yyaxis right
plot(t_cycle_extraUrb,vel_cycle_extraUrb,'Color', [0 0 0 0.3]);
set(ylabel('Velocit\`a [m/s]'),'Interpreter','latex');
yyaxis left
plot(t_cycle_extraUrb,SoC_opt_extraUrb,ColorMode="manual",Color="#A2142F");
hold on
plot(t_cycle_extraUrb,SoC_ECMS_extraUrb,ColorMode="manual",Color="#77AC30",LineStyle="-");
hold on
plot(t_cycle_extraUrb,SoC_CAECMS_extraUrb,ColorMode="manual",Color="#0072BD",LineStyle="-");
hold on
plot(t_cycle_extraUrb,SoC_DAECMS_extraUrb,ColorMode="manual",Color="#D95319",LineStyle="-");
hold on
yline(SoC_target,'--');
text(100,SoC_target+0.002,'SoC target')
set(ylabel('SoC'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
set(legend('PMP','ECMS','Continuous A-ECMS','Discrete A-ECMS',Location='bestoutside'),'Interpreter','latex');
title('Ciclo extra-urbano','FontWeight','bold');

subplot(3,1,3)
yyaxis right
plot(t_cycle_comb,vel_cycle_comb,'Color', [0 0 0 0.3]);
set(ylabel('Velocit\`a [m/s]'),'Interpreter','latex');
yyaxis left
plot(t_cycle_comb,SoC_opt_comb,ColorMode="manual",Color="#A2142F");
hold on
plot(t_cycle_comb,SoC_ECMS_comb,ColorMode="manual",Color="#77AC30",LineStyle="-");
hold on
plot(t_cycle_comb,SoC_CAECMS_comb,ColorMode="manual",Color="#0072BD",LineStyle="-");
hold on
plot(t_cycle_comb,SoC_DAECMS_comb,ColorMode="manual",Color="#D95319",LineStyle="-");
hold on
yline(SoC_target,'--');
text(100,SoC_target+0.002,'SoC target')
set(ylabel('SoC'),'Interpreter','latex');
set(xlabel('Time/sec'),'Interpreter','latex');
set(legend('PMP','ECMS','Continuous A-ECMS','Discrete A-ECMS',Location='bestoutside'),'Interpreter','latex');
title('Ciclo combinato','FontWeight','bold');
set(gcf,'Position',[200 100 1200 700]);
%%


% Consumo carburante solo ICE, PMP, ECMS, CA-ECMS, DA-ECMS
figure
b = bar(X_axis,Y_axis_1,'grouped','FaceColor','flat');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')

xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = string(b(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom')

xtips4 = b(4).XEndPoints;
ytips4 = b(4).YEndPoints;
labels4 = string(b(4).YData);
text(xtips4,ytips4,labels4,'HorizontalAlignment','center','VerticalAlignment','bottom')

xtips5 = b(5).XEndPoints;
ytips5 = b(5).YEndPoints;
labels5 = string(b(5).YData);
text(xtips5,ytips5,labels5,'HorizontalAlignment','center','VerticalAlignment','bottom')
set(gca,'XTick',1:3,'XTickLabel',{'Urbano','Extra-urbano', 'Combinato'},'TickLabelInterpreter','latex');
set(gca,'YTick',5,'YTickLabel',{'Consumo carburante [lt]'},'YTickLabelRotation',90,'TickLabelInterpreter','latex');
ylim([0 9]);
set(legend('Solo ICE','PMP','ECMS','Continuous A-ECMS','Discrete A-ECMS',Location='bestoutside'),'Interpreter','latex');
set(gcf,'Position',[200 100 1200 500]);
%%

% Plot risparmio carburante al traguardo PMP, ECMS, CA-ECMS, DA-ECMS
figure
b2 = bar(X_axis,Y_axis_2,'grouped','FaceColor','flat');
xtips1 = b2(1).XEndPoints;
ytips1 = b2(1).YEndPoints;
labels1 = string(b2(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')

xtips2 = b2(2).XEndPoints;
ytips2 = b2(2).YEndPoints;
labels2 = string(b2(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')

xtips3 = b2(3).XEndPoints;
ytips3 = b2(3).YEndPoints;
labels3 = string(b2(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom')

xtips4 = b2(4).XEndPoints;
ytips4 = b2(4).YEndPoints;
labels4 = string(b2(4).YData);
text(xtips4,ytips4,labels4,'HorizontalAlignment','center','VerticalAlignment','bottom')
set(gca,'XTick',1:3,'XTickLabel',{'Urbano','Extra-urbano', 'Combinato'},'TickLabelInterpreter','latex');
set(gca,'YTick',20,'YTickLabel',{'Risparmio carburante [\%]'},...
    'YTickLabelRotation',90,'TickLabelInterpreter','latex');
ylim([0 35]);
set(legend('PMP','ECMS','Continuous A-ECMS','Discrete A-ECMS',Location='bestoutside'),'Interpreter','latex');
set(gcf,'Position',[200 100 1200 500]);
