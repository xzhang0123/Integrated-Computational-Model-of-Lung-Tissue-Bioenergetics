
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
exp_dLAC0=[16 32.7 53.4 75 95.4]; 
exp_dLAC05=[16.8 34.3 54 78.7 100.1];
exp_dLAC2=[7.41 15.8 26.9 37 49.3];
exp_dLAC5=[4.77 12.5 19.9 27.7 35.8];
exp_dLAC20=[3.11 9.54 16.6 24.4 31.5];
%--------------------------------------
global  Tem Flow F_con R_con Vr Vb Vc Vm Vi...
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
   bufferpH=7.4;
    %------------
    Flow   =   12e-3;     %mL/min  
Vr  =   55e-3;      %mL
Vb  =   0.66e-3; %mL
Vcell  =  0.67e-3; %mL  

Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}   Faraday 's constant 
Tem=310.15; %K      37 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 
closed_system=0;
Para=ones(1,51);


  %%  Define t_step and t_final
t_step      =   0.1;   %min
%% Run Simulation
time0=20; 
time1=100;
tic
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);

LAC_STEP=[0 .5e-3 2e-3 5e-3 20e-3 ];
PYR_STEP=10*0.1*LAC_STEP;
GLC_Mass=zeros(5,time1/t_step+1);
LAC_Mass=zeros(5,time1/t_step+1);
for i=1:1:length(LAC_STEP)
    IC=Set_Initial_Concentrations;
IC(iGLCr)=10e-3;
IC(iLACr)=LAC_STEP(i); 
[T00,C00] = ode15s(@odeq,[0:0.1:time0],IC,options,Para);
IC1=C00(end,:);
[T,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
GLC_Mass(i,:)=1e6*Vr*(C(1,iGLCr)-C(:,iGLCr))/Lung_dw; % unit:umol/g dry weight
LAC_Mass(i,:)=1e6*Vr*(C(:,iLACr)-C(1,iLACr))/Lung_dw; %

end
set(figure(1),'Units','inches','Position',[0.2 0.1 5 4]) 
figure(1)

h2f1=plot(T_exp,exp_dLAC0,'^m','LineWidth',line_width,'MarkerSize',marker_size,'MarkerFaceColor','m')
hold on
h3f1=plot( T_exp,exp_dLAC05,'<k','LineWidth',line_width,'MarkerSize',marker_size,'MarkerFaceColor','k')

h4f1=plot(T_exp,exp_dLAC2,'ob','LineWidth',line_width,'MarkerSize',marker_size,'MarkerFaceColor','b')

h5f1=plot(T_exp,exp_dLAC5,'sr','LineWidth',line_width,'MarkerSize',marker_size,'MarkerFaceColor','r')
h6f1=plot(T_exp,exp_dLAC20,'*g','LineWidth',line_width,'MarkerSize',marker_size,'MarkerFaceColor','g')

h1f1=plot(T,GLC_Mass(1,:),'m',T,GLC_Mass(2,:),'k',T,GLC_Mass(3,:),'b',T,GLC_Mass(4,:),'r',T,GLC_Mass(5,:),'g',...
    'LineWidth',line_width,'MarkerSize',marker_size)
hold off
title('D.  GLC consumed')
box off
legend([h2f1 h3f1 h4f1 h5f1 h6f1],'0.0 mM LACr','0.5 mM LACr','2.0 mM LACr','5.0 mM LACr','20  mM LACr')
legend boxoff
xlabel('Recirculation time (min)')
ylabel('\mumol/g dry weight')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)