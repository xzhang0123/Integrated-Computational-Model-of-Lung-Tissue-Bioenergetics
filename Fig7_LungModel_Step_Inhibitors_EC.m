
clear all
close all
clc
format long
%--------------
size_index=dlmread('text_size.txt');
text_size=size_index(1);
text_size2=size_index(2);
line_width=size_index(3);
marker_size=size_index(4);
%%  Parameter Setup

global  Tem Flow F_con R_con Vr Vb Vc Vm Vi ...
   iGLCr iPYRr iLACr iPir iHr...
   iGLCb iPYRb iLACb iPib iHb...
   iGLCc iG6Pc iF6Pc iF16BPc iGAPc iBPGc iPEPc iPYRc iLACc iPG6c...
   iR5Pc iMALc iOXAc iCITc iaKGc iSUCc iFUMc iGLUc iASPc iPic...
   iAMPc iADPc iATPc iNADHc iNADc iNADPHc iNADPc iGSSG iGSH iH2O2...
    iHc   iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iGLUm...
    iASPm iNADm iNADHm iACOAm iCOAm iUQm iUQH2m  iPim iADPm iATPm...
    iFADm iFADH2m iHm iCytCoxi iCytCred iHi idPsim idPsip iO2 iR123e iR123m
    %------------------------------------------------
        %reservior
    iGLCr=1; iPYRr=2; iLACr=3 ;iPir=4; iHr=5;
    %blood 
    iGLCb=6; iPYRb=7; iLACb=8;iPib=9; iHb=10;
    %cytosol
    iGLCc=11 ;iG6Pc=12; iF6Pc=13; iF16BPc=14; iGAPc=15 ;iBPGc=16; iPEPc=17; iPYRc=18; 
    iLACc=19; iPG6c=20 ;iR5Pc=21; iMALc=22; iOXAc=23;
    iCITc=24; iaKGc=25; iSUCc=26; iFUMc=27; iGLUc=28; iASPc=29; iPic=30;
    iAMPc=31; iADPc=32; iATPc=33; iNADHc=34; iNADc=35;
    iNADPHc=36; iNADPc=37; iGSSG=38; iGSH=39; iH2O2=40; iHc=41;
    iPYRm=42; iOXAm=43; iCITm=44; iaKGm=45; iSCAm=46; iSUCm=47; iFUMm=48;
    iMALm=49; iGLUm=50; iASPm=51; iNADm=52; iNADHm=53; iACOAm=54; iCOAm=55;
    iUQm=56; iUQH2m=57;  iPim=58; iADPm=59; iATPm=60; iFADm=61; iFADH2m=62; iHm=63;
    %inter-membrane
    iCytCoxi=64; iCytCred=65; iHi=66;
    %other
    idPsim=67; idPsip=68; iO2=69; iR123e=70; iR123m=71;
    %------------
    Flow   =   12e-3;     %mL/min  
Lung_dw=0.227;
Vr  =   30e-3;      %mL
Lung_volume=1.6e-3;    %1.6mL 
Vb  =   0.66e-3; %mL
Vcell  =  0.67e-3; %mL  
UQss=1;
Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
 Tem=310.15; %K      30 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
IC=Set_Initial_Concentrations;
IC(iGLCc)=5.5e-3;
Para=ones(1,51);
  %%  Define t_step and t_final
t_step      =   0.1;   %min
%% Run Simulation
time0=10; %Length of State 2  (time before adding ADP)
time1=10;
time2=60;
time3=10;
tic
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);
Inhibitor_factor=[ 1 :-0.1:0.1    0.06   0.05 0.04 0.03 0.02 0.01 0];  %
%Inhibitor_factor=[ 1 :-0.1:0.1    0.06   0.05 0.04 0.03 ];  %
%  leak=1:2:2*length(Inhibitor_factor);
for i=1:1:length(Inhibitor_factor)
    
[T00,C00] = ode15s(@odeq,[0:t_step:time0],IC,options,Para);
IC1=C00(end,:);
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
Para1=Para;
Para1(27)=Inhibitor_factor(i)*Para(27);%CI 
%Para1(29)=Inhibitor_factor(i)*Para(29);%CIII
%Para1(45)=leak(i)*Para(45);%CI 
IC2=C1(end,:);
[T2,C2] = ode15s(@odeq,[0:t_step:time2],IC2,options,Para1);
%------calculate fluxes
Tfluxes=zeros(length(C2(:,1)),16)';
Rfluxes=zeros(length(C2(:,1)),31)';
for istep=1:1:(length(C2(1:end,1)))
    RTfluxes(:,istep)=fluxes(C2(istep,:),Para1);
end
Rfluxes=RTfluxes(1:31,:);
Tfluxes=RTfluxes(32:47,:);
CV(i)=1e6*(Rfluxes(31,end))/Lung_dw; 
ANT(i)=1e6*(Tfluxes(11,end))/Lung_dw; 
R_PPP(i)=1e6*abs(Rfluxes(9,end))/Lung_dw; %unit: umol/min/g dry weight 
T_LAC(i)=1e6*abs(Tfluxes(3,end))/Lung_dw; %unit: umol/min/g dry weight 
T_PYR(i)=1e6*abs(Tfluxes(2,end))/Lung_dw; %unit: umol/min/g dry weight 
T_PYRcm(i)=1e6*abs(Tfluxes(9,end))/Lung_dw; %unit: umol/min/g dry weight 
T_GLC(i)=1e6*abs(Tfluxes(1,end))/Lung_dw; %unit: umol/min/g dry weight 
LAC_PYR_ratio(i)=T_LAC(i)./T_PYR(i);
OCR(i)=1e6*Rfluxes(30,end)/Lung_dw;
% CI(i)=Rfluxes
NADHc_ratio(i)=C2(end,iNADHc)./C2(end,iNADc);
NADHm_ratio(i)=C2(end,iNADHm)./C2(end,iNADm);
T_NADH(i)=Tfluxes(15,end);
%calculate ATP ADP amount and Energy charge 
Mass_ATP=C2(:,iATPm)*Vm+C2(:,iATPc)*Vc;
Mass_ADP=C2(:,iADPm)*Vm+C2(:,iADPc)*Vc;
Mass_AMP=C2(:,iAMPc)*Vc;

MATP(i)=Mass_ATP(end);
MADP(i)=Mass_ADP(end);
MAMP(i)=Mass_AMP(end);
Mtotal(i)=MATP(i)+MADP(i)+MAMP(i);
EC(i)=(Mass_ATP(end)+0.5*Mass_ADP(end))./(Mass_ATP(end)+Mass_ADP(end)+Mass_AMP(end));
ALA(i)=1e6*abs(Tfluxes(16,end))/Lung_dw;
ddPsim(i)=C2(end,idPsim);
ATPc(i)=C2(end,iATPc);
ADPc(i)=C2(end,iADPc);
AMPc(i)=C2(end,iAMPc);
Pic(i)=C2(end,iPic);
PPTL(i)=ATPc(i)/(ADPc(i)*Pic(i));
ATPm(i)=C2(end,iATPm);
ADPm(i)=C2(end,iADPm);
Pim(i)=C2(end,iPim);
end

EC_ANP=dlmread('EC_ANP.txt');
ANP_DATA=EC_ANP(:,1);
EC_DATA=EC_ANP(:,2);
XEC_AMP=ANP_DATA(1:21);
XEC_ADP=ANP_DATA(22:47);
XEC_ATP=ANP_DATA(48:71);
DATAY_AMP=EC_DATA(1:21);
DATAY_ADP=EC_DATA(22:47);
DATAY_ATP=EC_DATA(48:71);
set(figure(1),'Units','inches','Position',[0.2 0.1 5 4]) 
set(figure(2),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(1)
h1f1=plot(XEC_ATP,DATAY_ATP,'r*','LineWidth',line_width,'markersize',marker_size)

hold on
h2f1=plot(XEC_ADP,DATAY_ADP,'g^','LineWidth',line_width,'markersize',marker_size)
h3f1=plot(XEC_AMP,DATAY_AMP,'bd','LineWidth',line_width,'markersize',marker_size)
h4f1=plot(EC,MATP./Mtotal,'r',EC,MADP./Mtotal,'g',EC,MAMP./Mtotal,'b','LineWidth',line_width,'markersize',marker_size)
hold off
legend([h1f1 h2f1 h3f1 ],'ATP','ADP','AMP')
box off
xlim([0.1 1])
legend boxoff
xlabel('Energy charge')
ylabel('Adenine nucleotide fraction')
  set(gcf,'color','w')
set(gca,'Fontsize',16)
EC_LAC=dlmread('EC_LAC.txt');
EXP_LAC=EC_LAC(:,1)/15;% change unit from umol/15 min/g dw to umol/min/g dw
EXP_EC=EC_LAC(:,2);
figure(2)
plot(EC,T_LAC,'r',EXP_EC,EXP_LAC,'r^','LineWidth',line_width,'markersize',marker_size)
xticks([0.2:0.2:1])
%title('LAC Production VS Energy Charge')
xlabel({'Energy charge'})
legend('Model','Data')
legend boxoff
ylabel({'LAC production rate';'(\mumol/min/g dry weight)'})
box off
  set(gcf,'color','w')
set(gca,'Fontsize',16)
ylim([0 3.5])
xlim([0.1 1])
set(figure(4),'Units','inches','Position',[3 0.1 5 4]) 
figure(4)
plot(Inhibitor_factor,EC,'r','LineWidth',line_width,'markersize',marker_size)
%title('LACPYR ratio and NADH ratio')
xlabel('Fraction of normal CI activity')
ylabel('Energy charge')
box off
  set(gcf,'color','w')
set(gca,'Fontsize',16)
ylim([0 1])
% 
% set(figure(5),'Units','inches','Position',[0.2 0.1 5 4]) 
% figure(5)
% plot(T_LAC/1,(36.03+R_con*Tem*log(PPTL)),'r','LineWidth',line_width,'markersize',marker_size)
% %title('LAC Production VS \DeltaG_{ATPase}')
% xlabel({'LAC production rate';'(\mumol/min/g dry weight)'})
% ylabel('-\DeltaG_{ATPase} (KJ/mol)')
% xlim([0 3.5])
% % legend('Model','Data')
% % legend boxoff
% box off
%   set(gcf,'color','w')
% set(gca,'Fontsize',16)
% xlim([0 3.5])
% ylim([0 1])
% 
% set(figure(6),'Units','inches','Position',[0.2 0.1 5 4]) 
% figure(6)
% plot(EC,PPTL,'r','LineWidth',line_width,'markersize',marker_size)
% title('Phosphate potential VS Energy Charge')
% xlabel('Energy charge')
% ylabel('Phosphate potential ')
% box off
%   set(gcf,'color','w')
% set(gca,'Fontsize',16)
% % 
% set(figure(7),'Units','inches','Position',[0.2 0.1 5 4]) 
% figure(7)
% 
% h71=plot(Inhibitor_factor,CV,Inhibitor_factor,ANT,'LineWidth',line_width,'markersize',marker_size)
% %title('Energy Charge VS \DeltaG_{ATPase}')
% xticks([0.0:0.2:1])
% box off 
% 
% xlabel('Fraction of normal CIII activity')
% ylabel('reaction flux (\mumol/min/g dry weight))')
% legend('CV','ANT')
% legend boxoff
%   set(gcf,'color','w')
% set(gca,'Fontsize',16)

set(figure(8),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(8)

yyaxis left
h81=plot(Inhibitor_factor,-(-36.03-R_con*Tem*log(PPTL)),'LineWidth',line_width,'markersize',marker_size)
%title('Energy Charge VS \DeltaG_{ATPase}')
xticks([0.0:0.2:1])
xlabel('Fraction of normal CI activity')
ylabel('-\DeltaG_{ATPase} (KJ/mol)')
ylim([44 56])
hold on
yyaxis right
h82=plot(Inhibitor_factor,EC,'LineWidth',line_width,'markersize',marker_size)
ylabel('Energy charge')
box off
  set(gcf,'color','w')
set(gca,'Fontsize',16)
ylim([0 1])

set(figure(9),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(9)

plot(Inhibitor_factor,1e3*Pim,'k',Inhibitor_factor,1e3*ATPm,'r',Inhibitor_factor,1e3*ADPm,'g','LineWidth',line_width)
box off
ylabel('Concentrations (mM)')
xlabel('Fraction of normal CI activity')
legend('Pi_m','ATP_m','ADP_m')
legend boxoff
  set(gcf,'color','w')
set(gca,'Fontsize',16)



set(figure(10),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(10)
yyaxis left
h101=plot(Inhibitor_factor,OCR,'LineWidth',line_width)
xticks([0.0:0.2:1])
ylabel('OCR (\mumol/min/g dry weight)')

yyaxis right
h102=plot(Inhibitor_factor,T_LAC,'LineWidth',line_width)
box off
xlabel ('Fraction of normal CI activity')
ylabel('LAC (\mumol/min/g dry weight)')
  set(gcf,'color','w')
set(gca,'Fontsize',16)



set(figure(11),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(11)
plot(Inhibitor_factor,ddPsim)
% plot(100*Inhibitor_factor,-(-36.03-R_con*Tem*log(PPTL)))
% title('OCR VS -\DeltaG_{ATPase}')
% ylabel ('\DeltaG_{ATPase}')
% xlabel('OCR')
%   set(gcf,'color','w')
% set(gca,'Fontsize',16)



set(figure(12),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(12)
h8=plot(Inhibitor_factor,1e3*Pic,'k',Inhibitor_factor,1e3*ATPc,'r',Inhibitor_factor,1e3*ADPc,'g',Inhibitor_factor,1e3*AMPc,'b','LineWidth',line_width)
ylabel('Concentrations (mM)')
xlabel('Fraction of normal CI activity')
xlim([0.0 1])
ylim([0 3.5])
% legend(h8,{'Pi','ATP','ADP'})
legend([ h8],{'Pi_c','ATP_c','ADP_c','AMP_c'})
legend boxoff
ylabel('Concentrations (mM)')
box off
  set(gcf,'color','w')
set(gca,'Fontsize',16)


% 
% 
% set(figure(13),'Units','inches','Position',[0.2 0.1 5 4]) 
% figure(13)
% yyaxis left
% h131=plot(Inhibitor_factor,1e3*ATPc,'-r',Inhibitor_factor,1e3*ADPc,'-g',Inhibitor_factor,1e3*AMPc,'-b','LineWidth',line_width)
% ylabel('ANP Concentrations (mM)')
% xlabel('Fraction of normal CI activity')
% xticks([0.0:0.2:1])
% xlim([0.0 1])
% ylim([0 3.5])
% % legend(h8,{'Pi','ATP','ADP'})
% 
% yyaxis right
% h132=plot(Inhibitor_factor,1e3*Pic,'k','LineWidth',line_width)
% legend boxoff
% ylabel('Pi Concentration (mM)')
% legend([h131; h132],{'ATP_c','ADP_c','AMP_c','Pi_c'})
% legend boxoff
% box off
%   set(gcf,'color','w')
% set(gca,'Fontsize',16)
% 
% 
% 
% set(figure(14),'Units','inches','Position',[0.2 0.1 5 4]) 
% figure(14)
% yyaxis left
% h141=plot(Inhibitor_factor,1e3*ATPm,'-r',Inhibitor_factor,1e3*ADPm,'-g','LineWidth',line_width)
% ylabel('ANP Concentrations (mM)')
% xlabel('Fraction of normal CI activity')
% xticks([0.0:0.2:1])
% xlim([0.0 1])
% % legend(h8,{'Pi','ATP','ADP'})
% yyaxis right
% h142=plot(Inhibitor_factor,1e3*Pim,'k','LineWidth',line_width)
% legend boxoff
% ylabel('Pi Concentration (mM)')
% legend([h141; h142],{'ATP_m','ADP_m','Pi_m'})
% legend boxoff
% box off
%   set(gcf,'color','w')
% set(gca,'Fontsize',16)