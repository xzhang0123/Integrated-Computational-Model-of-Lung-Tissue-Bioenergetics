
clear all
close all
clc
format long
%--------------
size_index=dlmread('text_size.txt');
text_size=1.2*size_index(1);
text_size2=size_index(2);
line_width=size_index(3);
marker_size=0.7*size_index(4);
%%  Parameter Setup
inhibitor=input('Choose an inhibitor:\n  1: CI inhibitor ROT;\n  3: CIII inhibitor AA;\n  4: CIV inhibitor CO or KCN;\n  6: Uncoupler;\n  7: PGI glycolysis inhibitor;\n  8: MA shuttle inhibitor\n  9: ROT+AA CI and CIII inhibitors\n')

%--------
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
  

Vr  =   55e-3;      %mL
Lung_volume=1.6e-3;    %1.6mL 
Vb  =   0.66e-3; %mL
Vcell  =  0.67e-3; %mL  

Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
Tem=310.15; %K      30 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 
ROTi=1;%CI inhibited by ROT : percentage
closed_system=0;
IC=Set_Initial_Concentrations;
 Para=ones(1,51);

  %%  Define t_step and t_final
t_step      =   0.05;   %min
%% Run Simulation
time0=20; %Length of State 2  (time before adding ADP)
time1=2.8;
time2=1.5;
time3=10.5;
tic
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);


[T00,C00] = ode15s(@odeq,[0:t_step:time0],IC,options,Para);
IC1=C00(end,:);
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
 Para1=Para;

if inhibitor==1
   Para1(27)=0.1*0.15*Para(27);%CI inhibitor
elseif inhibitor==3
     Para1(29)=0.0001*Para(29);%CIII inhibitor   AA
elseif inhibitor==4 
    Para1(30)=0.0001*Para(30);% KCN, CIV inhibitor
elseif inhibitor==6    
Para1(45)=5*Para(45);%Uncoupler
elseif inhibitor==7 
Para1(2)=0.00001*Para(2);%PGI inhibitor
elseif inhibitor==8
Para1(46)=0.05*Para(46);%MA shuttle inhibitor
elseif inhibitor==9
Para1(27)=0.03*Para(27);%CI inhibitor
Para1(29)=0.0*Para(29);%CIII inhibitor   AA

end
IC2=C1(end,:);
[T2,C2] = ode15s(@odeq,[0:t_step:time2],IC2,options,Para1);
IC3=C2(end,:);
Para2=Para1;
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
HEX=1e9*Rfluxes(1,:);
PGI=1e9*Rfluxes(2,:);
PFK=1e9*Rfluxes(3,:);
G3PF=1e9*Rfluxes(4,:);
G3PD=1e9*Rfluxes(5,:);
PHK=1e9*Rfluxes(6,:);
PK=1e9*Rfluxes(7,:);
LD=1e9*Rfluxes(8,:);
PPP1=1e9*Rfluxes(9,:);
PPP2=1e9*Rfluxes(10,:);
GSH1=1e9*Rfluxes(11,:);
GSH2=1e9*Rfluxes(12,:);
ATPase=1e9*Rfluxes(13,:);
AK=1e9*Rfluxes(14,:);
MDH2=1e9*Rfluxes(15,:);
GOT2=1e9*Rfluxes(16,:);
PDH=1e9*Rfluxes(17,:);
CITS=1e9*Rfluxes(18,:);
CITDH=1e9*Rfluxes(19,:);
AKGDH=1e9*Rfluxes(20,:);
SCAS=1e9*Rfluxes(21,:);
NDK=1e9*Rfluxes(22,:);
SUCDH=1e9*Rfluxes(23,:);
FH=1e9*Rfluxes(24,:);
MDH1=1e9*Rfluxes(25,:);
GOT1=1e9*Rfluxes(26,:);
CI=1e9*Rfluxes(27,:);
CII=1e9*Rfluxes(28,:);
CIII=1e9*Rfluxes(29,:);
CIV=1e9*Rfluxes(30,:);
CV=1e9*Rfluxes(31,:);

Tr1=1e9*Tfluxes(1,:);
Tr2=1e9*Tfluxes(2,:);
Tr3=1e9*abs(Tfluxes(3,:));
Tr4=1e9*Tfluxes(4,:);
Tr5=1e9*Tfluxes(5,:);
Tr6=1e9*Tfluxes(6,:);
Tr7=1e9*Tfluxes(7,:);
Tr8=1e9*Tfluxes(8,:);
Tr9=1e9*Tfluxes(9,:);
Tr10=1e9*Tfluxes(10,:);
Tr11=1e9*Tfluxes(11,:);
Tr12=1e9*Tfluxes(12,:);
Tr13=1e9*Tfluxes(13,:);
Tr14=1e9*Tfluxes(14,:);
Tr15=1e9*Tfluxes(15,:);
Tr16=1e9*Tfluxes(16,:);
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


MNADH=(Vm*C(:,iNADHm)+Vc*C(:,iNADHc))/(Vm+Vc);
MFAD=Vm*C(:,iFADm)/Vm;

%--------- backgroud noise
if inhibitor==1
bg_NADH_rot=1.6*3.12e-5;
MNADH=MNADH+bg_NADH_rot;
elseif inhibitor==4
   bg_NADH_KCN=2.4e-5; 
   MNADH=MNADH+bg_NADH_KCN;
   MFAD=MFAD+200*bg_NADH_KCN;
   elseif inhibitor==6
   %bg_NADH_UCP=0.4e-5; 
    bg_NADH_UCP=1.2e-5; 
   MNADH=MNADH+bg_NADH_UCP;
   MFAD=MFAD+60*bg_NADH_UCP;
end

%---------------------------
MNADH=MNADH/MNADH(1);
MFAD=MFAD/MFAD(1);
%-------------------------------------------------------
EXP_FADNADH=xlsread('NADH_FAD_DATA.xlsx',1);
EXP_time_ROT=EXP_FADNADH(1:49,1)/60;
EXP_ROT_NADH=EXP_FADNADH(1:49,2);
EXP_ROT_FAD=EXP_FADNADH(1:49,3);
EXP_ROT_RATIO=EXP_ROT_NADH./EXP_ROT_FAD;

EXP_time_KCN=EXP_FADNADH(1:49,5)/60;
EXP_KCN_NADH=EXP_FADNADH(1:49,6);
EXP_KCN_FAD=EXP_FADNADH(1:49,7);
EXP_KCN_RATIO=EXP_KCN_NADH./EXP_KCN_FAD;

EXP_time_PCP=EXP_FADNADH(1:49,9)/60;
EXP_PCP_NADH=EXP_FADNADH(1:49,10);
EXP_PCP_FAD=EXP_FADNADH(1:49,11);
EXP_PCP_RATIO=EXP_PCP_NADH./EXP_PCP_FAD;
set(figure(7),'Units','inches','Position',[0.2 0.1 5 4]) 
if inhibitor==1
figure(7)
plot(EXP_time_ROT,EXP_ROT_NADH,'b*',T,MNADH,'k','LineWidth',2,'markersize',marker_size)
box off
title('ROT effect')
legend('Data','Model')
legend boxoff
xlabel('Recirculation time (min)','FontName','Times New Roman')
ylabel('NADH')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlim([0 6])
ylim([0.7 1.3])
elseif inhibitor==4
figure(7)
plot(EXP_time_KCN,EXP_KCN_NADH,'b*',T,MNADH,'k','LineWidth',2,'markersize',marker_size)
title('KCN effect')
box off
% xlabel('Recirculation time (min)','FontName','Times New Roman')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlim([0 6])
ylim([0.7 1.3])
xlabel('Recirculation time (min)','FontName','Times New Roman')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
elseif inhibitor==6
figure(7)
plot(EXP_time_PCP,EXP_PCP_NADH,'b*',T,MNADH,'k','LineWidth',2,'markersize',marker_size)
title('PCP effect')
box off
% legend('NADH (Data)','NADH (Model)')
% legend boxoff
xlabel('Recirculation time (min)','FontName','Times New Roman')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlim([0 6])
ylim([0.7 1.3])
end