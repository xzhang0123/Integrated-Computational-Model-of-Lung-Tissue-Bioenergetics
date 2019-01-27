
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
Lung_dw=0.227; %each lung has 0.227g dry weight (1.33/5.87)
T_exp=20:20:100;
T_exp2=[19 39 59 79 99];
exp_glu05=[2.36 5.35 8.55 12.6 15.4]; %Unit: 
lbexp_glu05=[1 1 2 3 2];
ubexp_glu05=[1 1 2 2 2];

exp_glu56=[10.8 23.5 35.7 56.2 69.5];
lbexp_glu56=[2 5 3.5 2 6];
ubexp_glu56=[2 5 3.5 6 10];

exp_glu15=[4.67 11 19.5 27.7 36.4];
lbexp_glu15=[2.09 1.4 0.8 1.7 2.8];
ubexp_glu15=[0.81 1.8 2.1 1.1 3.6];

exp_glu3=[7.2 16.3 28.8 41.9 53.5];
lbexp_glu3=[2.06 0.5 3.3 6 9];
ubexp_glu3=[2.45 2.5 3.2 6 9];

exp_glu10=[13.5 24.5 46.2 61.4 83.2];
lbexp_glu10=[2 5 3 2 1];
ubexp_glu10=[2 5 3 8 7];
% %---------------lactate data------------
exp_lac05_exo=[13.1 24.8 30.4 32.6 33]*0.18;  %Subtract endogenase part
exp_lac15_exo=[20.5 32.6 44.7 52.9 61]*0.36;
exp_lac3_exo=[20.6 33.5 45.7 57.2 68.1]*0.41;
exp_lac56_exo=[20 38.3 54.2 70.7 84.8]*0.46;
exp_lac10_exo=[22.6 48.8 69 88.7 111]*0.75;

% %--------
%---------------lactate data------------
exp_lac05=[13.1 24.8 30.4 32.6 33];  %Subtract endogenase part
lbexp_lac05=[1.5 1.5 2 3 3];
ubexp_lac05=[1.5 1.8 2 3 3];

exp_lac3=[20.6 33.5 45.7 57.2 68.1];
lbexp_lac3=[2 3 2 5 2];
ubexp_lac3=[2 3 2 5 2];

exp_lac15=[20.5 32.6 44.7 52.9 61];
lbexp_lac15=[2 3 6.4 6.4 5.2];
ubexp_lac15=[2 3 6 5 6.5];

exp_lac56=[20 38.3 54.2 70.7 84.8];
lbexp_lac56=[2 3 2 5 5.2];
ubexp_lac56=[2 3 5 5 5.2];

exp_lac10=[22.6 48.8 69 88.7 111];
lbexp_lac10=[2 5 8 10 13.1];
ubexp_lac10=[2 5 9 10 13.1];
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
   substrates=0;
   bufferpH=7.4;
    %------------
    Flow   =   12e-3;     %mL/min  


Vr  =   55e-3;      %mL

Vb  =   0.66e-3; %mL
Vcell  =  0.67e-3; %mL  
Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
  Tem=310.15; %K      37 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]

Para=ones(1,51);


  %%  Define t_step and t_final
t_step      =   1;   %min
%% Run Simulation
time0=20; 
time1=100;
tic
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);
GLC_STEP=[0.5e-3 1.5e-3 3e-3 5.6e-3 10e-3 ];
GLC_Mass=zeros(5,time1/t_step+1);
LAC_Mass=zeros(5,time1/t_step+1);

for i=1:1:length(GLC_STEP)
tic
IC=Set_Initial_Concentrations;
    IC(iGLCr)=GLC_STEP(i);
    IC(iGLCb)=GLC_STEP(i);
[T00,C00] = ode15s(@odeq,[0:0.1:time0],IC,options,Para);
IC1=C00(end,:);
% Para(29)=0.5*Para(29);%CIII
[T,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
toc
GLC_Mass(i,:)=1e6*Vr*(C(1,iGLCr)-C(:,iGLCr))/Lung_dw; %change unit from mol/lung to umol/g dw
LAC_Mass(i,:)=1e6*Vr*(C(:,iLACr)-C(1,iLACr))/Lung_dw; %change unit from mol to umol
%--------------------------------------
Tfluxes=zeros(length(C(:,1)),16)';
Rfluxes=zeros(length(C(:,1)),31)';
for istep=1:1:(length(C(1:end,1)))
    RTfluxes(:,istep)=fluxes(C(istep,:),Para);
end
Rfluxes=RTfluxes(1:31,:);
Tfluxes=RTfluxes(32:47,:);
%----------------------------------
SIM_LAC_total(i)=abs(1e6*Tfluxes(3,end))/Lung_dw;
SIM_LAC_exo(i)=Rfluxes(7,end)./(Tfluxes(16,end)+Rfluxes(7,end))*SIM_LAC_total(i); 
SIM_TOTAL_ALA(i)=Tfluxes(16,end)/Lung_dw;
SIM_LAC_endo(i)=abs(1e6*Tfluxes(3,end))/Lung_dw-SIM_LAC_exo(i);
end

set(figure(1),'Units','inches','Position',[0.2 0.1 5 4]) 
set(figure(2),'Units','inches','Position',[0.2 0.1 5 4]) 

figure(1)
errorbar(T_exp,exp_glu05,lbexp_glu05,ubexp_glu05,'gs', 'MarkerFaceColor','g','MarkerEdgeColor','g', 'MarkerSize',marker_size, 'LineWidth',line_width)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
box off
hold on
errorbar(T_exp,exp_glu15,lbexp_glu15,ubexp_glu15,'rs', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',marker_size, 'LineWidth',line_width)
errorbar(T_exp,exp_glu3,lbexp_glu3,ubexp_glu3,'sb', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerSize',marker_size, 'LineWidth',line_width)
errorbar(T_exp,exp_glu56,lbexp_glu56,ubexp_glu56,'ks', 'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',marker_size, 'LineWidth',line_width)
errorbar(T_exp2,exp_glu10,lbexp_glu10,ubexp_glu10,'ms', 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',marker_size, 'LineWidth',line_width)
plot(T,GLC_Mass(1,:),'g',T,GLC_Mass(2,:),'r',T,GLC_Mass(3,:),'b',T,GLC_Mass(4,:),'k',T,GLC_Mass(5,:),'m','LineWidth',line_width)

title('A.  GLC consumed')
legend('0.5 mM GLCr','1.5 mM GLCr','3.0 mM GLCr','5.6 mM GLCr','10 mM  GLCr')
legend boxoff
xlabel('Recirculation time (min)')
ylabel('\mumol/g dry weight')
hold off


figure(2)

errorbar(T_exp,exp_lac05,lbexp_lac05,ubexp_lac05,'dg', 'MarkerFaceColor','g','MarkerEdgeColor','g', 'MarkerSize',marker_size, 'LineWidth',line_width)
hold on
ylim([0 150])
errorbar(T_exp,exp_lac15,lbexp_lac15,ubexp_lac15,'dr', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerSize',marker_size, 'LineWidth',line_width)
errorbar(T_exp2,exp_lac3,lbexp_lac3,ubexp_lac3,'db', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerSize',marker_size, 'LineWidth',line_width)
errorbar(T_exp,exp_lac56,lbexp_lac56,ubexp_lac56,'sk', 'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',marker_size, 'LineWidth',line_width)
errorbar(T_exp2,exp_lac10,lbexp_lac10,ubexp_lac10,'^m', 'MarkerFaceColor','m','MarkerEdgeColor','m', 'MarkerSize',marker_size, 'LineWidth',line_width)
plot(T,LAC_Mass(1,:),'g',T,LAC_Mass(2,:),'r',T,LAC_Mass(3,:),'b',T,LAC_Mass(4,:),'k',T,LAC_Mass(5,:),'m','LineWidth',line_width)
box off
title('B.  Total LAC produced')
legend('0.5 mM GLCr','1.5 mM GLCr','3.0 mM GLCr','5.6 mM GLCr','10 mM  GLCr')
legend boxoff
xlabel('Recirculation time (min)')
ylabel('\mumol/g dry weight')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
hold off


GLC_C=[  1.5 5.6 10];
LAC_endo=[ 43.4 39.6 31.9 ]/100;
LAC_endo_SE=[ 2.6 7 4.1]/100;
LAC_total=[ 68 83.1 111 ]/100;
LAC_total_SE=[ 3.5 5.2 3.1]/100
LAC_exo=LAC_total-LAC_endo;
LAC_exo_SE=[ 3.5 5.2 3.1]/100

set(figure(3),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(3)
h71_exo=plot(GLC_STEP*1e3,SIM_LAC_exo,'b','linewidth',2)
box off
hold on
h72_exo= errorbar(GLC_C,LAC_exo,LAC_exo_SE,'d','marker','square','markersize',marker_size,...
                 'markeredgecolor','b', 'MarkerFaceColor','b',...
                'color','b','linewidth',2)

h71_endo=plot(GLC_STEP*1e3,SIM_LAC_endo,'r','linewidth',2)
h72_endo= errorbar(GLC_C,LAC_endo,LAC_endo_SE,'d','marker','^','markersize',marker_size,...
                 'markeredgecolor','r', 'MarkerFaceColor','r',...
                'color','r','linewidth',2)
            %----------          
h71_total=plot(GLC_STEP*1e3,SIM_LAC_total,'k','linewidth',2)
hold on
h72_total= errorbar(GLC_C,LAC_total,LAC_total_SE,'d','marker','o','markersize',marker_size,...
                 'markeredgecolor','k', 'MarkerFaceColor','k',...
                'color','k','linewidth',2)
hold off
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
xlabel('Glucose concentration (mM)')
ylabel('\mumol/min/g dry weight')
title('C.  LAC production')
       legend([h72_exo h72_endo h72_total],'Exogenous LAC','Endogenous LAC','Total LAC')
       legend boxoff
   hold off   
  ylim([0 1.5])
   