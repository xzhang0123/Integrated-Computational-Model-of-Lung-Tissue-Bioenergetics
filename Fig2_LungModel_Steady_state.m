
clear all
close all
clc
format long
%------define text size for figures------------
text_size=1.2*16;
text_size2=18;
line_width=2;
marker_size=12;
%%  Parameter Setup
%--------

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
 Lung_dw=0.227; 
    Flow   =   12e-3;     %mL/min  

Vr  =   50e-3;      %mL
Vb  =   0.66e-3; %0.66mL   blood volume
Vcell  =  0.67e-3; %0.67mL  

Vm=1/51.0714*Vcell;%  calculate m c i volumes
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}    Faraday 's constant
 Tem=310.15; %K      37 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 
IC=Set_Initial_Concentrations;
 IC(iGLCr)=5.5e-3;
Para=ones(1,51);  %Adjustable parameters


  %%  Define t_step and t_final
t_step      =   0.1;   %min
%% Run Simulation
time0=20; %Length of
time1=100;
tic
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);


[T00,C00] = ode15s(@odeq,[0:0.1:time0],IC,options,Para);
IC1=C00(end,:);
[T,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
Tfluxes=zeros(length(C(:,1)),16)';
Rfluxes=zeros(length(C(:,1)),31)';
for istep=1:1:(length(C(1:end,1)))
    RTfluxes(:,istep)=fluxes(C(istep,:),Para);
end
Rfluxes=RTfluxes(1:31,:)./Lung_dw;
Tfluxes=RTfluxes(32:47,:)./Lung_dw;
HEX=1e6*Rfluxes(1,:);
PGI=1e6*Rfluxes(2,:);
PFK=1e6*Rfluxes(3,:);
G3PF=1e6*Rfluxes(4,:);
G3PD=1e6*Rfluxes(5,:);
PHK=1e6*Rfluxes(6,:);
PK=1e6*Rfluxes(7,:);
LD=1e6*Rfluxes(8,:);
PPP1=1e6*Rfluxes(9,:);
PPP2=1e6*Rfluxes(10,:);
GSR=1e6*Rfluxes(11,:);
GP=1e6*Rfluxes(12,:);
ATPase=1e6*Rfluxes(13,:);
AK=1e6*Rfluxes(14,:);
MDH2=1e6*Rfluxes(15,:);
GOT2=1e6*Rfluxes(16,:);
PDH=1e6*Rfluxes(17,:);
CITS=1e6*Rfluxes(18,:);
CITDH=1e6*Rfluxes(19,:);
AKGDH=1e6*Rfluxes(20,:);
SCAS=1e6*Rfluxes(21,:);
NDK=1e6*Rfluxes(22,:);
SUCDH=1e6*Rfluxes(23,:);
FH=1e6*Rfluxes(24,:);
MDH1=1e6*Rfluxes(25,:);
GOT1=1e6*Rfluxes(26,:);
CI=1e6*Rfluxes(27,:);
CII=1e6*Rfluxes(28,:);
CIII=1e6*Rfluxes(29,:);
CIV=1e6*Rfluxes(30,:);
CV=1e6*Rfluxes(31,:);
GLUT=1e6*Tfluxes(1,:);
PYRT=1e6*Tfluxes(2,:);
LACT=1e6*Tfluxes(3,:);
Tr4=1e6*Tfluxes(4,:);
Tr5=1e6*Tfluxes(5,:);
Tr6=1e6*Tfluxes(6,:);
Tr7=1e6*Tfluxes(7,:);
Tr8=1e6*Tfluxes(8,:);
Tr9=1e6*Tfluxes(9,:);
Tr10=1e6*Tfluxes(10,:);
Tr11=1e6*Tfluxes(11,:);
Tr12=1e6*Tfluxes(12,:);
Tr13=1e6*Tfluxes(13,:);
Tr14=1e6*Tfluxes(14,:);
Tr15=1e6*Tfluxes(15,:);
ALA=1e6*Tfluxes(16,:);
ATPm_rate=1e6*Rfluxes(21,:)+1e6*Rfluxes(31,:);
ATPc_rate=-1e6*Rfluxes(1,:)-1e6*Rfluxes(3,:)+1e6*Rfluxes(6,:)+1e6*Rfluxes(7,:);
CO2_rate=1e6*Rfluxes(17,:)+1e6*Rfluxes(19,:)+1e6*Rfluxes(20,:);  %unit is umol/g dw/min

MATP=Vm*C(:,iATPm)+Vc*C(:,iATPc);
MADP=Vm*C(:,iADPm)+Vc*C(:,iADPc);
MAMP=Vc*C(:,iAMPc);
EC=(MATP+0.5*MADP)./(MATP+MADP+MAMP);

GLYProfile=[1e6*Rfluxes(1:8,end); 1e6*Tfluxes(15,end)];
PPPProfile=1e6*Rfluxes([9 11:12],end);
TCAProfile=1e6*Rfluxes([17:21 23 25],end);
ETCProfile=1e6*Rfluxes(27:31,end);
GLY_P_Profile=1e6*Rfluxes([1:3 5 7:9],end);

Tprofilebc=1e6*abs(Tfluxes(1:3,end));
% set(figure(5),'Units','inches','Position',[0.2 0.1 6 4]) 
% set(figure(6),'Units','inches','Position',[0.2 0.1 6 4])
% set(figure(7),'Units','inches','Position',[0.2 0.1 6 4]) 
% set(figure(8),'Units','inches','Position',[0.2 0.1 6 4])
set(figure(9),'Units','inches','Position',[0.2 0.1 6 4]) 

figure(5)
bar(GLYProfile,'b')
set(gca,'XTickLabel',{'HEX','PGI','PFK','GAPF','GAPDH','PHK','PK','LDH','MAS'})
   set(gcf,'color','w')
set(gca,'Fontsize',0.75*text_size)
ylabel('Fluxes (nmol/min)','Fontsize',text_size)
title('Glycolysis reaction fluxes','Fontsize',text_size)
box off
PPP_index=[1];
GLYP_index=[7];
mean_PPPrate=[0.158*0.115/Lung_dw  ];  %
se_PPPrate=[0.002/Lung_dw];
set(figure(14),'Units','inches','Position',[0.2 0.1 6 4]) 
figure(14)
h1_GLYP=bar(GLY_P_Profile,'b');
set(gca,'XTickLabel',{'HEX','PGI','PFK','GAPDH','PK','LDH','PPP'})
   set(gcf,'color','w')
set(gca,'Fontsize',0.75*text_size)
ylabel('Fluxes (nmol/min)','Fontsize',text_size)
title('Glycolysis and pentose phosphate cycle','Fontsize',text_size)
box off
hold on
h2_GLYP=errorbar(GLYP_index,mean_PPPrate,se_PPPrate,'d','marker','*','markersize',5,...
                 'markeredgecolor','r',...
                'color','r','linewidth',3,'CapSize',10);
            legend([h2_GLYP h1_GLYP],'Data','Simulation')
         legend boxoff
         hold off

%total oxygen consumption rate=2500*0.227=568, CIV=2100*0.227=488
% mean_ETCrate=[2*488 363*1000/60*0.227];  %oxygen consumption rate=568, CIV activity=2*OCR
ETC_index=[ 4 ];%ATPchange unit from umol/h/g dwt to nmol/min/lung
mean_ETCrate=[2*0.488 ]/Lung_dw;  %oxygen consumption rate=568, CIV activity=2*OCR
se_ETCrate=[2*0.050 ]/Lung_dw;
%--------------
figure(8)
h1_ETC=bar(ETCProfile,'b');
set(gca,'XTickLabel',{'CI','CII','CIII','CIV','CV'})
box off
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)
ylabel('Fluxes (\mumol/min/g dry weight)','Fontsize',text_size)
title('ETC reaction fluxes','Fontsize',text_size)
hold on
h2_ETC=errorbar(ETC_index,mean_ETCrate,se_ETCrate,'d','marker','*','markersize',10,...
                 'markeredgecolor','r',...
                'color','r','linewidth',3,'CapSize',10)
            legend([h2_ETC h1_ETC],'Data','Simulation')
         legend boxoff
       hold off
 %------------------------------------
TGLY_index=[1 2 3 ];
mean_Trate=[158 192.5/10 192.5 ]/1e3/Lung_dw;
se_Trate=[15.8,19.25/10,19.25]/1e3/Lung_dw;
figure(9)
h1_T=bar(Tprofilebc,0.85,'b');
set(gca,'XTickLabel',{'GLC','PYR','LAC'})
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)
box off
ylabel('Fluxes (\mumol/min/g dry weight)','Fontsize',text_size)
title('A.  Transport fluxes','Fontsize',text_size)
ylim([0 1.2])
hold on
h2_T=errorbar(TGLY_index,mean_Trate,se_Trate,'d','marker','*','markersize',10,...
                 'markeredgecolor','r',...
                'color','r','linewidth',3,'CapSize',10);
             legend([h2_T h1_T],'Data','Model')
             legend boxoff
                set(gcf,'color','w')
 hold off
 energy_index=[1 2 3];
 mean_energy=[5.66/1.17 5.66/0.31 (5.66+0.5*0.31)/(5.66+1.17+0.31)];
 se_energy=0.15*mean_energy;
 Energy_profile=[MATP(end)./MADP(end) MATP(end)./MAMP(end) EC(end)];
       

 set(figure(14),'Units','inches','Position',[0.2 0.1 3 4]) 
figure(14)
h1_PPP2=bar(1,PPPProfile(1),1.2,'b');
set(gca,'XTickLabel',{'G6PDH'})
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)
ylabel('Fluxes (\mumol/min/g dry weight)','Fontsize',text_size)
title('B.  PPP','Fontsize',text_size)
box off
ylim([0 0.1])
hold on
h2_PPP2=errorbar(PPP_index(1),mean_PPPrate(1),se_PPPrate(1),'d','marker','*','markersize',10,...
                 'markeredgecolor','r',...
                'color','r','linewidth',3,'CapSize',10);
            legend([h2_PPP2 h1_PPP2],'Data','Model')
         legend boxoff
         hold off
         
         
 set(figure(15),'Units','inches','Position',[0.2 0.1 3 4]) 
figure(15)
h1_OCR=bar(1,0.5*ETCProfile(4),0.6,'b');
set(gca,'XTickLabel',{'OCR'})
box off
   set(gcf,'color','w')
set(gca,'Fontsize',text_size)
ylabel('Fluxes (\mumol/min/g dry weight)','Fontsize',text_size)
title('C.  OCR','Fontsize',text_size)
hold on
h2_OCR=errorbar(1,0.5*mean_ETCrate,0.5*se_ETCrate,'d','marker','*','markersize',10,...
                 'markeredgecolor','r',...
                'color','r','linewidth',3,'CapSize',10);
            legend([h2_ETC h1_OCR],'Data','Simulation')
         legend boxoff
       hold off
       
  
 set(figure(16),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(16)
ATPP=[ATPc_rate(end) ATPm_rate(end)];
ATPP_data=[1.03 6.05];
h1_ATPP=bar([1 2],ATPP,0.6,'b');
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
hold on
box off
title('D.  ATP production','Fontsize',text_size)
ylabel('Fluxes (\mumol/min/g dry weight)','Fontsize',text_size)
h2_ATPP=errorbar([1 2],ATPP_data,0.0*ATPP_data,'d','marker','*','markersize',10,...
                 'markeredgecolor','r',...
                'color','r','linewidth',3,'CapSize',10);
            legend([h2_ATPP h1_ATPP],'Data','Model')
         legend boxoff
         
 
       hold off
       set(gca,'XTickLabel',{'Cytosolic','Mitochondrial'})
       
      