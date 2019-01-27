
clear all
close all
clc
format long
%--------------
size_index=dlmread('text_size.txt');
text_size=1*size_index(1);
text_size2=size_index(2);
line_width=size_index(3);
marker_size=0.8*size_index(4);
%%  Parameter Setup

EXP_Inh=xlsread('Inhibitor_Fisher_Chance_1976');
inhibitor=input('Choose a inhibitor:\n  1: CI inhibitor ROT;\n  3: CIII inhibitor AA;\n  4: CIV inhibitor CO or KCN;\n  6: Uncoupler;\n  7: PGI glycolysis inhibitor;\n  8: MA shuttle inhibitor\n  9: Hyperoxia\n')
Lung_dw=0.227;
Time_deglu=EXP_Inh(1:3,1);
deglu_LAC=EXP_Inh(1:3,2);
SE_deglu_LAC=EXP_Inh(1:3,4);
deglu_PYR=EXP_Inh(1:3,3);
SE_deglu_PYR=EXP_Inh(1:3,5);
%-----------------------------
Time_CO=EXP_Inh(7:10,1);
CO_LAC=EXP_Inh(7:10,2)/1e3;  %change unit from nmol to umol
SE_CO_LAC=EXP_Inh(7:10,4)/1e3;
CO_PYR=EXP_Inh(7:10,3)/1e3;
SE_CO_PYR=EXP_Inh(7:10,5)/1e3;
%------------------
Time_AOA=EXP_Inh(13:20,1)-10.6;
AOA_LAC=EXP_Inh(13:20,2);
SE_AOA_LAC=EXP_Inh(13:20,4);
AOA_PYR=EXP_Inh(13:20,3);
SE_AOA_PYR=EXP_Inh(13:20,5);
global  Tem Flow F_con R_con Vr Vb Vc Vm Vi ROTi closed_system...
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
    %LUNG WEIGH 1.6 mL,  Vb is 10%  Vc is 90%, Vmito=~2% cells Vi=~10% of Mito

Vr  =   30e-3;      %mL
%Lung_volume=1.6e-3;    %1.6mL 
Vb  =   0.66e-3; %0.66 mL
Vcell  =  0.67e-3; %mL  

Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
 Tem=310.15; %K      37 oC
 %Tem=298.15; %K      25 oC   25 degree
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 
 % Ve  =   2*1e-3; %L 


ROTi=1;%CI inhibited by ROT : percentage
closed_system=0;
IC=Set_Initial_Concentrations;
 IC(iGLCr)=5.5e-3;
%  if inhibitor==4
%   IC(iO2)=95/21*IC(iO2);
%  end
% IC(iPYRc)=5e-3;
% IC(iMALc)=5e-3;
%IC(iNADHc)=3e-3;
Para=ones(1,51);
  %%  Define t_step and t_final
t_step      =  0.1;   %min
%% Run Simulation
time0=100; 
time1=20;
time2=10;
time3=15;
tic
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);
%options = odeset('NonNegative',[1:71]);

[T00,C00] = ode15s(@odeq,[0:t_step:time0],IC,options,Para);
IC1=C00(end,:);

[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
Para1=Para;

%--------------------------------------------
if inhibitor==1
    Para1(27)=0.15*Para(27);%CI inhibitor
    %ROT non-competative with ROT, reduces Vmax
elseif inhibitor==3
     Para1(29)=0.0001*Para(29);%CIII inhibitor   AA
elseif inhibitor==4 
%Para1(30)=0.01*Para(30);% CO, CIV inhibitor
 % Para1(30)=1*0.003*Para(30);% CO, CIV inhibitor
   Para1(51)=1e7*Para(51);% CO, CIV inhibitor
     % Para1(51)=1e6*Para(51);% CO, CIV inhibitor
elseif inhibitor==6    
Para1(45)=5*Para(45);%Uncoupler
elseif inhibitor==7 
Para1(2)=0.0000*Para(2);%PGI inhibitor GLY
elseif inhibitor==8
Para1(46)=0.005*Para(46);%MA shuttle inhibitor
elseif inhibitor==9
    %Hyperoxia
Para1(27)=0.23*Para(27);%CI inhibitor
Para1(28)=0.36*Para(28);%CII inhibitor   
Para1(23)=0.36*Para(23);%SDH
elseif inhibitor==10
    %Hyperoxia
Flow=0;
Para1(49)=0*Para(49);
%Para1(30)=0.1*Para(30);% CO, CIV inhibitor
end
%----------------------------------
IC2=C1(end,:);
% IC2(iGLCr)=5.5e-3;
% IC2(iPYRr)=0;
% IC2(iLACr)=0;
%IC2(iO2)=0.01*IC2(iO2);
[T2,C2] = ode15s(@odeq,[0:t_step:time2],IC2,options,Para1);
IC3=C2(end,:);
%Para2=Para;   %delete inhibitor
Para2=Para1;   %keep inhibitor
if inhibitor==8
Para2(46)=0.02*Para(46);%MA shuttle inhibitor
end
%Para2=Para1; %keep inhibitor
[T3,C3] = ode15s(@odeq,[0:t_step:time3],IC3,options,Para2);
T=[T1; T2(2:end)+time1; T3(2:end)+time1+time2;];

C=[C1; C2(2:end,:);C3(2:end,:);];
Tfluxes=zeros(length(C(:,1)),16)';
Rfluxes=zeros(length(C(:,1)),31)';
for istep=1:1:(length(C1(1:end,1)))
    RTfluxes1(:,istep)=fluxes(C1(istep,:),Para);
end
for istep=1:1:(length(C2(1:end,1)))
    RTfluxes2(:,istep)=fluxes(C2(istep,:),Para1);
end
for istep=1:1:(length(C3(1:end,1)))
    RTfluxes3(:,istep)=fluxes(C3(istep,:),Para2);
end
RTfluxes=[RTfluxes1 RTfluxes2(:,2:end) RTfluxes3(:,2:end)];
Rfluxes=RTfluxes(1:31,:);
Tfluxes=RTfluxes(32:47,:);
clear t t_step t_final
HEX=1e6*Rfluxes(1,:)/Lung_dw;
PGI=1e6*Rfluxes(2,:)/Lung_dw;
PFK=1e6*Rfluxes(3,:)/Lung_dw;
G3PF=1e6*Rfluxes(4,:)/Lung_dw;
G3PD=1e6*Rfluxes(5,:)/Lung_dw;
PHK=1e6*Rfluxes(6,:)/Lung_dw;
PK=1e6*Rfluxes(7,:)/Lung_dw;
LD=1e6*Rfluxes(8,:)/Lung_dw;
PPP1=1e6*Rfluxes(9,:)/Lung_dw;
PPP2=1e6*Rfluxes(10,:)/Lung_dw;
GSH1=1e6*Rfluxes(11,:)/Lung_dw;
GSH2=1e6*Rfluxes(12,:)/Lung_dw;
ATPase=1e6*Rfluxes(13,:)/Lung_dw;
AK=1e6*Rfluxes(14,:)/Lung_dw;
MDH2=1e6*Rfluxes(15,:)/Lung_dw;
GOT2=1e6*Rfluxes(16,:)/Lung_dw;
PDH=1e6*Rfluxes(17,:)/Lung_dw;
CITS=1e6*Rfluxes(18,:)/Lung_dw;
CITDH=1e6*Rfluxes(19,:)/Lung_dw;
AKGDH=1e6*Rfluxes(20,:)/Lung_dw;
SCAS=1e6*Rfluxes(21,:)/Lung_dw;
NDK=1e6*Rfluxes(22,:)/Lung_dw;
SUCDH=1e6*Rfluxes(23,:)/Lung_dw;
FH=1e6*Rfluxes(24,:)/Lung_dw;
MDH1=1e6*Rfluxes(25,:)/Lung_dw;
GOT1=1e6*Rfluxes(26,:)/Lung_dw;
CI=1e6*Rfluxes(27,:)/Lung_dw;
CII=1e6*Rfluxes(28,:)/Lung_dw;
CIII=1e6*Rfluxes(29,:)/Lung_dw;
CIV=1e6*Rfluxes(30,:)/Lung_dw;
CV=1e6*Rfluxes(31,:)/Lung_dw;

Tr1=1e6*Tfluxes(1,:)/Lung_dw;
Tr2=1e6*abs(Tfluxes(2,:))/Lung_dw; %b to c is the primary direction
Tr3=1e6*abs(Tfluxes(3,:))/Lung_dw;
Tr4=1e6*Tfluxes(4,:)/Lung_dw;
Tr5=1e6*Tfluxes(5,:)/Lung_dw;
Tr6=1e6*Tfluxes(6,:)/Lung_dw;
Tr7=1e6*Tfluxes(7,:)/Lung_dw;
Tr8=1e6*Tfluxes(8,:)/Lung_dw;
Tr9=1e6*Tfluxes(9,:)/Lung_dw;
Tr10=1e6*Tfluxes(10,:)/Lung_dw;
Tr11=1e6*Tfluxes(11,:)/Lung_dw;
Tr12=1e6*Tfluxes(12,:)/Lung_dw;
Tr13=1e6*Tfluxes(13,:)/Lung_dw;
Tr14=1e6*Tfluxes(14,:)/Lung_dw;
Tr15=1e6*Tfluxes(15,:)/Lung_dw;
Tr16=1e6*Tfluxes(16,:)/Lung_dw;

ATPm_rate=1e6*Rfluxes(21,:)+1e6*Rfluxes(31,:);
ATPc_rate=-1e6*Rfluxes(1,:)-1e6*Rfluxes(3,:)+1e6*Rfluxes(6,:)+1e6*Rfluxes(7,:);
%--------MASS BALANCE
 ATPM = (C(:,iAMPc)+C(:,iATPc)+C(:,iADPc))*Vc+(C(:,iATPm)+C(:,iADPm))*Vm;
% % CoA/ACoA/SCA
 CoAM = (C(:,iCOAm)+C(:,iACOAm)+C(:,iSCAm))*Vm;
% % TCA compounds CIT AKG SCA SUC MAL OAA  ASP GLU FUM
TCAM = (C(:,iOXAc)+C(:,iMALc)+C(:,iCITc)+C(:,iaKGc)+C(:,iSUCc)+C(:,iFUMc)+C(:,iGLUc)+C(:,iASPc))*Vc+(C(:,iGLUm)+C(:,iASPm)+C(:,iOXAm)+C(:,iCITm)+C(:,iaKGm)+C(:,iSCAm)+C(:,iSUCm)+C(:,iFUMm)+C(:,iMALm))*Vm;

mass_balance.adenosine      = ATPM;
mass_balance.coenzymeA      = CoAM;
mass_balance.TCAcompounds   = TCAM;

fprintf('Variance in mass of adenosine:\t\t%f micromol, %f percent of initial total\n',1e6*max(mass_balance.adenosine)-1e6*min(mass_balance.adenosine),(max(mass_balance.adenosine)-min(mass_balance.adenosine))/mass_balance.adenosine(1)*100)
fprintf('Variance in mass of CoA:\t\t%f micromol, %f percent of initial total\n',1e6*max(mass_balance.coenzymeA)-1e6*min(mass_balance.coenzymeA),(max(mass_balance.coenzymeA)-min(mass_balance.coenzymeA))/mass_balance.coenzymeA(1)*100)
fprintf('Variance in mass of TCA compounds:\t%f micromol, %f percent of initial total\n',1e6*max(mass_balance.TCAcompounds)-1e6*min(mass_balance.TCAcompounds),(max(mass_balance.TCAcompounds)-min(mass_balance.TCAcompounds))/mass_balance.TCAcompounds(1)*100)

%-----------------
set(figure(1),'Units','inches','Position',[0.2 0.1 12 10]) 
set(figure(2),'Units','inches','Position',[0.2 0.1 12 10]) 
set(figure(3),'Units','inches','Position',[0.2 0.1 12 10]) 
set(figure(4),'Units','inches','Position',[0.2 0.1 12 10]) 
figure(1)
subplot(4,4,1)
plot(T,C(:,iGLCr),T,C(:,iPYRr),T,C(:,iLACr),T,C(:,iPir),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Reservior Concentrations')
legend('GLCr','PYRr','LACr','Pir')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,2)
plot(T,C(:,iGLCb),T,C(:,iPYRb),T,C(:,iLACb),T,C(:,iPib),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Concentrations in blood region')
legend('GLCb','PYRb','LACb','Pib')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,3)
plot(T,C(:,iGLCc),T,100*C(:,iPYRc),T,C(:,iLACc),T,C(:,iPic),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('GLCc','100*PYRc','LACc','Pic')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,4)
plot(T,C(:,iG6Pc),T,C(:,iF6Pc),T,C(:,iF16BPc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('G6Pc','F6Pc','F18BPc')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,5)
plot(T,C(:,iGAPc),T,C(:,iBPGc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('GAPc','BPGc')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,6)
plot(T,C(:,iPEPc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('PEP')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,7)
plot(T,C(:,iNADHc),T,C(:,iNADc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('NADHc','NADc')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,8)
plot(T,C(:,iATPc),T,C(:,iADPc),T,C(:,iAMPc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('ATPc','ADPc','AMPc')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,9)
plot(T,C(:,iMALc),T,C(:,iSUCc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('MALc','SUCc')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,10)
plot(T,C(:,iASPc),T,C(:,iOXAc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('ASPc','OXAc')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)

subplot(4,4,11)
plot(T,C(:,iGLUc),T,C(:,iaKGc),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('GLUc','aKGc')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)

subplot(4,4,12)
plot(T,C(:,iPYRc),T,C(:,iPYRm),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('cytosol concentrations')
legend('PYRc','PYRm')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)



%-----------------------------------------
figure(2)
subplot(4,4,1)
plot(T,C(:,iPYRm),T,C(:,iCITm),T,C(:,iSUCm),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('PYRm','CITm','SUCm')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,2)
plot(T,C(:,iNADHm),T,C(:,iNADm),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('NADH','NAD')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,3)
plot(T,C(:,iFADH2m),T,C(:,iFADm),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('FADH2','FAD')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,4)
plot(T,C(:,idPsim),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('dPsi')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,5)
plot(T,C(:,iADPm),T,C(:,iATPm),T,C(:,iPim),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('ADPm','ATPm','Pim')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,6)
plot(T,C(:,iCOAm),T,C(:,iACOAm),T,C(:,iSCAm),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('COA','ACOA','SCA')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,7)
plot(T,C(:,iMALm),T,C(:,iSUCm),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('MALm','SUCm')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,8)
plot(T,C(:,iGLUm),T,C(:,iASPm),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('GLUm','ASPm')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,9)
plot(T,C(:,iCytCoxi),T,C(:,iCytCred),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('CytCo','CytCr')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,10)
plot(T,C(:,iUQm),T,C(:,iUQH2m),'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito concentrations')
legend('UQm','UQH2m')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,11)
Mass_ATP=C(:,iATPm)*Vm+C(:,iATPc)*Vc;
Mass_ADP=C(:,iADPm)*Vm+C(:,iADPc)*Vc;
Mass_AMP=C(:,iAMPc)*Vc;
EC=(Mass_ATP+0.5*Mass_ADP)./(Mass_ATP+Mass_ADP+Mass_AMP);
plot(T,EC,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Energy Charge')
legend('Energy Charge')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
%
subplot(4,4,12)
plot(T,Mass_ATP,T,Mass_ADP,T,Mass_AMP,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Energy Charge')
legend('ATP','ADP','AMP')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,13)
plot(T,Mass_ATP./Mass_ADP,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('ATP/ADP ratio')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,14)
plot(T,HEX,T,PPP1,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('HEX and PPP')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
legend('HEX','PPP')
      
subplot(4,4,15)
 %convert unit from nmol/g dw/min to umol/h/g dry wt
 ATPc_rate= ATPc_rate*60;
  ATPm_rate= ATPm_rate*60;
 plot(T,ATPc_rate,T,ATPm_rate,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('ATP production rate')
legend('ATPc','ATPm')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)


figure(3)
subplot(4,4,1)
plot(T,HEX,T,PGI,T,PFK,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Mito reaction fluxes (nmol/min)')
legend('HEX','PGI','PFK')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)

subplot(4,4,2)
plot(T,G3PF,T,G3PD,T,PHK,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('G3PF','G3PD','PHK')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)

subplot(4,4,3)
plot(T,PK,T,LD,T,PPP1,T,PPP2,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('PK','LD','PPP1','PPP2')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,4)
plot(T,ATPase,T,AK,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('ATPase','AK')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,5)
plot(T,MDH2,T,GOT2,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('MD2','GOT2')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)

subplot(4,4,6)
plot(T,PDH,T,CITS,T,CITDH,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('PDH','CITS','CITDH')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)


subplot(4,4,7)
plot(T,AKGDH,T,SCAS,T,SUCDH,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('AKGDH','SCAS','SUCDH')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)


subplot(4,4,8)
plot(T,AKGDH,T,SCAS,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('AKGDH','SCAS')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,9)
plot(T,CI,T,CII,T,CIII,T,CIV,T,CV,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('reaction fluxes (nmol/min)')
legend('CI','CII','CIII','CIV','CV')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
ylabel('Flux (nmol/min)','FontName','Times New Roman', 'FontSize',10)


figure(4)
subplot(4,4,1)
plot(T,Tr1,T,Tr2,T,Tr3,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Transport fluxes (nmol/min)')
legend('GLCT b-c','PYRT b-c','LAC b-c')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,2)
plot(T,Tr5,T,Tr6,T,Tr7,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Transport fluxes (nmol/min)')
legend('SUC-PI','MAL-PI','MAL-AKG')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)

subplot(4,4,3)
plot(T,Tr8,T,Tr9,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Transport fluxes (nmol/min)')
legend('T(TCC)','T(PYRH cm)')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,4)
plot(T,Tr10,T,Tr11,T,Tr12,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Transport fluxes (nmol/min)')
legend('T(PIC)','T(ANT)','T(GLUH)')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,5)
plot(T,Tr13,T,Tr14,T,Tr15,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Transport fluxes (nmol/min)')
legend('T(GAE)','T(LEAK)','T(NADH)')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)
subplot(4,4,6)
plot(T,Tr16,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Transport fluxes (nmol/min)')
legend('T16')
xlabel('Time (min)','FontName','TiTes New Roman', 'FontSize',10)

subplot(4,4,7)
plot(T,Tr4,'LineWidth',1.2)
   set(gcf,'color','w')
set(gca,'Fontsize',10)
title('Transport fluxes (nmol/min)')
legend('PI(b-c)')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',10)



if inhibitor== 7
set(figure(7),'Units','inches','Position',[2 1 5 7]) 
figure(7)
subplot(2,1,1)
plot(T,Tr2,'r','LineWidth',line_width,'markersize',marker_size)
title('Effect of deoxyglucose (PGI inhibitor)')
hold on
h2=errorbar(Time_deglu,deglu_PYR,SE_deglu_PYR,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',marker_size, 'LineWidth',line_width)
legend('PYR (Model)','PYR (Normalized Data)')
legend boxoff
xlabel('Time (min)')
ylabel('nmol/min/Lung')
box off
   set(gcf,'color','w')
set(gca,'Fontsize',16)
xlim([0 30])
hold off
subplot(2,1,2)
plot(T,Tr3,'b','LineWidth',line_width,'markersize',marker_size)
hold on 
h2=errorbar(Time_deglu,deglu_LAC,SE_deglu_LAC,'sb', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerSize',marker_size, 'LineWidth',line_width)
xlim([0 30])
legend('LAC (Model)','LAC (Normalized Data)')
legend boxoff
xlabel('Time (min)')
ylabel('nmol/min/Lung')
box off
   set(gcf,'color','w')
set(gca,'Fontsize',16)
hold off
end
%% 
    if inhibitor ==4
% 
  set(figure(104),'Units','inches','Position',[3 3 5 4]) 
 figure(104)
yyaxis left
h2=errorbar(Time_CO-5,CO_PYR,SE_CO_PYR,'db', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerSize',marker_size, 'LineWidth',line_width)
%title({'C.  LAC and PYR Production';'(\mumol/min/g dry weight)'})
hold on
box off
h1=plot(T-15,Tr2,'b','LineWidth',line_width)
ylabel({'PYR production';'(\mumol/min/g dry weight)'})
xlim([0 20])
ylim([0 0.5])
set(gca,'YTick',0:0.1:0.5) 
set(gca,'XTick',0:5:20) 
yyaxis right
h4=errorbar(Time_CO-5,1.6*CO_LAC,SE_CO_LAC,'sr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',marker_size, 'LineWidth',line_width)
h3 =plot(T-15,Tr3,'r','LineWidth',line_width)
hold off
legend([ h4  h2],'LAC','PYR')
legend boxoff
xlabel('Recirculation time (min)')
ylabel({'LAC production';'(\mumol/min/g dry weight)'})
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)
ylim([0 3])
hold off
    end 
    %% 
if inhibitor==8
figure(8)
subplot(2,1,1)
plot(T,Tr2,'r','LineWidth',line_width)
xlim([0 30])
hold on
h2=errorbar(Time_AOA,AOA_PYR,SE_AOA_PYR,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',marker_size, 'LineWidth',line_width)
legend('PYR (Model)','PYR (Data)')
hold off
title('Effect of AOA')
box off
legend boxoff
xlabel('Time (min)')
ylabel('nmol/min/Lung')
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)

subplot(2,1,2)
plot(T,Tr3,'b','LineWidth',line_width)
xlim([0 30])
ylim([20 140])
hold on
h2=errorbar(Time_AOA,AOA_LAC,SE_AOA_LAC,'sb', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerSize',marker_size, 'LineWidth',line_width)
hold off
title('Effect of AOA')
legend('LAC (Model)','LAC (1.8*Data)')
box off
legend boxoff
xlabel('Time (min)')
ylabel('nmol/min/Lung')
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)
end
if inhibitor ==1 
    T_exp=[5 15];  %Normal condition and Inhibitor condition
    exp_LACPYR_ratio=[14 48];
    lbexp_LACPYR_ratio=[4 10];
    ubexp_LACPYR_ratio=[4 10];
      exp_LAC=[12.36 38.62]/15; %change unit from (umol/15min/g dw) to (umol/min/dw)
    lbexp_LAC=[1.64 3.14]/15;
    ubexp_LAC=[1.64 3.14]/15;
         exp_PYR=[0.88 0.87]/15;
    lbexp_PYR=[0.11 0.11]/15;
    ubexp_PYR=[0.14 0.14]/15;
    %-------------------------------
    exp_ATPADP_ratio=[5.66/1.17 2.34/1.83];
    lbexp_ATPADP_ratio=[0.10*5.66/1.17 0.10*2.34/1.83]
    ubexp_ATPADP_ratio=lbexp_ATPADP_ratio;
    %-----------------------------
     exp_ATPAMP_ratio=[5.66/0.31 2.34/1.73];
    lbexp_ATPAMP_ratio=[0.15*5.66/0.31 0.15*2.34/1.73]
    ubexp_ATPAMP_ratio=lbexp_ATPAMP_ratio;
    set(figure(101),'Units','inches','Position',[2 1 4 4]) 
      set(figure(102),'Units','inches','Position',[2 1 5 4]) 
        set(figure(103),'Units','inches','Position',[2 1 5 4]) 
    figure(101)
 
h1=plot(T-10,Mass_ATP./Mass_ADP,'r','LineWidth',line_width,'markersize',marker_size)
ylim([0 6])
xlim([0 20])
hold on 
h2=errorbar(T_exp,exp_ATPADP_ratio,lbexp_ATPADP_ratio,ubexp_ATPADP_ratio,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',marker_size, 'LineWidth',line_width)
hold off
box off
% legend([h1 h2],'Model prediction','Data')
% legend boxoff
ylabel('ATP/ADP ratio')
%title('A.  ATP/ADP ratio ')
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlabel('Recirculation time (min)')
figure(102)
yyaxis left
 h1=plot(T-10,Tr2,'b','LineWidth',line_width,'markersize',marker_size)
 set(gca,'YTick',0:0.1:0.5) 
 ylim([0 0.5])
 xlim([0 20])
 box off
 hold on 
 h3=errorbar(T_exp,exp_PYR,lbexp_PYR,ubexp_PYR,'sb', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerSize',marker_size, 'LineWidth',line_width)
ylabel({'PYR production';'(\mumol/min/g dry weight)'})
 yyaxis right
  ylim([0 3.5])
  h1=plot(T-10,Tr3,'r','LineWidth',line_width,'markersize',marker_size)
 h2=errorbar(T_exp,exp_LAC,lbexp_LAC,ubexp_LAC,'sr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',marker_size, 'LineWidth',line_width)
legend([h2,h3],'LAC','PYR')
legend boxoff
hold off
xlabel('Recirculation time (min)')
ylabel({'LAC production';'(\mumol/min/g dry weight)'})
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
%  figure(103)
%  h3=plot(T,Tr2,'r','LineWidth',line_width,'markersize',marker_size)
%  hold on
%  h2=errorbar(T_exp,exp_PYR,lbexp_PYR,ubexp_PYR,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',8, 'LineWidth',line_width)
% hold off
%     figure(109)
%     subplot(2,1,1)
%  h1=plot(T,Mass_ATP./Mass_AMP,'r','LineWidth',line_width,'markersize',marker_size)
%   title('ATP/AMP ratio (ROT)')
% hold on 
% h2=errorbar(T_exp,exp_ATPAMP_ratio,lbexp_ATPAMP_ratio,ubexp_ATPAMP_ratio,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',8, 'LineWidth',line_width)
% hold off
% legend([h1 h2],'Simulation','Data')
% legend boxoff
% box off
%    set(gcf,'color','w')
% set(gca,'Fontsize',16)
%   subplot(2,1,2)
%  plot(T,EC,'r','LineWidth',line_width,'markersize',marker_size)
%  box off
%  title('Energy Charge')
% xlabel('Time (min)')
% ylabel('EC')
%    set(gcf,'color','w')
% set(gca,'Fontsize',16)
end

figure(10)
 h1=plot(T,Tr3/Lung_dw,'b',T,Tr2/Lung_dw,'r','LineWidth',line_width,'markersize',marker_size)
  legend ('LAC','PYR')
  legend boxoff
  title('LAC and PYR production (ROT)')
xlabel('Time (min)')
ylabel('\mumol/min/g dry weight')
   set(gcf,'color','w')
set(gca,'Fontsize',16)
%-----------------------

if inhibitor ==3
    T_exp=[10 30];  %Normal condition and Inhibitor condition
    exp_LACPYR_ratio=[14 24.38];
    lbexp_LACPYR_ratio=[4 4];
    ubexp_LACPYR_ratio=[4 4];
    %-------------------------------
    exp_ATPADP_ratio=[5.66/1.17 1.68/1.58];
    lbexp_ATPADP_ratio=[0.15*5.66/1.17 0.15*1.68/1.58]
    ubexp_ATPADP_ratio=lbexp_ATPADP_ratio;
    %-----------------------------
     exp_ATPAMP_ratio=[5.66/0.31 1.68/2.84];
    lbexp_ATPAMP_ratio=[0.15*5.66/0.31 0.15*1.68/2.84]
    ubexp_ATPAMP_ratio=lbexp_ATPAMP_ratio;
    figure(8)
    subplot(2,1,1)
h1=plot(T,Mass_ATP./Mass_ADP,'r','LineWidth',line_width,'markersize',marker_size)
ylim([0 6])
hold on 
h2=errorbar(T_exp,exp_ATPADP_ratio,lbexp_ATPADP_ratio,ubexp_ATPADP_ratio,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',8, 'LineWidth',line_width)
hold off
box off
legend([h1 h2],'Simulation','Data')
legend boxoff
ylabel('ratio')
title('ATP/ADP ratio (AA)')
   set(gcf,'color','w')
set(gca,'Fontsize',16)
  subplot(2,1,2)
 h1=plot(T,Tr3./Tr2,'r','LineWidth',line_width,'markersize',marker_size)
  ylim([0 200])
 hold on 
h2=errorbar(T_exp,exp_LACPYR_ratio,lbexp_LACPYR_ratio,ubexp_LACPYR_ratio,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',8, 'LineWidth',line_width)
hold off
box off
% legend([h1 h2],'Simulation','Data')
% legend boxoff
 title('LAC/PYR ratio (AA)')

xlabel('Time (min)')
ylabel('ratio')
   set(gcf,'color','w')
set(gca,'Fontsize',16)
    figure(9)
    subplot(2,1,1)
 h1=plot(T,Mass_ATP./Mass_AMP,'r','LineWidth',line_width,'markersize',marker_size)
  title('ATP/AMP ratio (AA)')
hold on 
h2=errorbar(T_exp,exp_ATPAMP_ratio,lbexp_ATPAMP_ratio,ubexp_ATPAMP_ratio,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',8, 'LineWidth',line_width)
hold off
legend([h1 h2],'Simulation','Data')
legend boxoff
box off
   set(gcf,'color','w')
set(gca,'Fontsize',16)
  subplot(2,1,2)
 plot(T,EC,'r','LineWidth',line_width,'markersize',marker_size)
 box off
 title('Energy Charge (AA)')
xlabel('Time (min)')
ylabel('EC')
   set(gcf,'color','w')
set(gca,'Fontsize',16)
figure(10)
 h1=plot(T,Tr3/Lung_dw,'b',T,10*Tr2/Lung_dw,'r','LineWidth',line_width,'markersize',marker_size)
  legend ('LAC','10*PYR')
  legend boxoff
  title('LAC and PYR production (AA)')
xlabel('Time (min)')
ylabel('\mumol/min/g dry weight')
   set(gcf,'color','w')
set(gca,'Fontsize',16)
end