
clear all
close all
clc
format long
%--------------
size_index=dlmread('text_size.txt');
text_size=1.2*size_index(1);
text_size2=size_index(2);
line_width=size_index(3);
marker_size=size_index(4);
%%  Parameter Setup

global  Tem Flow F_con R_con Vr Vb Vc Vm Vi  ...
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
    Lung_dw=0.227;
    Flow   =   15e-3;     %mL/min  

Vr  =   40e-3;      %mL
Lung_volume=1.3e-3;    %1.3mL 
Vb  =   0.66e-3; %mL
Vcell  =  0.67e-3; %mL  

Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
 Tem=310.15; %K      30 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 


closed_system=0;
IC=Set_Initial_Concentrations;
Para=ones(1,51);

  %%  Define t_step and t_final
t_step      =   1;   %min
%% Run Simulation
time0=10; 
time1=10;
time2=80;
time3=10;
tic
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);
dPO2=[0.04:0.25:2.04]; % mmHg
dPO2=[dPO2 5 10 100 500 760];
%convert dPo2 to O2 concentration
dC_O2=0.003*dPO2/(22.4*1e6); %Unit: M
%  Para(27)=0.5*Para(27);%CI
for i=1:1:length(dPO2) 
[T00,C00] = ode15s(@odeq,[0:t_step:time0],IC,options,Para);
IC1=C00(end,:);
IC1(iO2)=dC_O2(i);
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
 
IC2=C1(end,:);
[T2,C2] = ode15s(@odeq,[0:t_step:time2],IC2,options,Para);
%------calculate fluxes
Tfluxes=zeros(length(C2(:,1)),16)';
Rfluxes=zeros(length(C2(:,1)),31)';
for istep=1:1:(length(C2(1:end,1)))
    RTfluxes(:,istep)=fluxes(C2(istep,:),Para);
end
Rfluxes=RTfluxes(1:31,:);
Tfluxes=RTfluxes(32:47,:);
R_PPP(i)=1e6*abs(Rfluxes(9,end))/Lung_dw; %unit: umol/min  absolute volue
T_LAC(i)=1e6*abs(Tfluxes(3,end))/Lung_dw; %unit: umol/min  absolute volue
T_PYR(i)=-1e6*(Tfluxes(2,end))/Lung_dw; %unit: umol/min  absolute volue
T_PYRcm(i)=1e6*abs(Tfluxes(9,end))/Lung_dw; %unit: umol/min  absolute volue
T_GLC(i)=1e6*abs(Tfluxes(1,end))/Lung_dw; %unit: umol/min  absolute volue
LAC_PYR_ratio(i)=T_LAC(i)./T_PYR(i);
NADH_ratio(i)=C2(end,iNADHc)./C2(end,iNADc);
%end of calculate fluxes-----
CLAC(i)=C2(end,iLACc);
CPYR(i)=C2(end,iPYRc);
Mass_ATP=C2(:,iATPm)*Vm+C2(:,iATPc)*Vc;
Mass_ADP=C2(:,iADPm)*Vm+C2(:,iADPc)*Vc;
Mass_AMP=C2(:,iAMPc)*Vc;

MATP(i)=Mass_ATP(end);
MADP(i)=Mass_ADP(end);
MAMP(i)=Mass_AMP(end);
Mtotal(i)=MATP(i)+MADP(i)+MAMP(i);
EC(i)=(Mass_ATP(end)+0.5*Mass_ADP(end))./(Mass_ATP(end)+Mass_ADP(end)+Mass_AMP(end));
Ratio_ATPADP(i)=MATP(i)./MADP(i);
end

FEXP_PO2= [0.95 0.2 0.05 0.01 0.001 0.00006 ]*(760);  
FEXP_LAC=[69.4 67.4 85.3 81.4 110 131 ]/60; %change unit umol/hr/g dry wt to umol/min/g dry wt
FEXP_LAC_SE=[5.3 3 4.6 4.7 7.9 6.6]/60;
FEXP_PYR=[7.9 7.4 7.9 6.7 6.6 6.6 ]/60;
FEXP_PYR_SE=[0.4 0.3 0.3 0.2 0.3 0.1]/60;
FEXP_RATIO_SE=[0.4 0.5 0.5 1.1 1.3 0.9 ]; 
FEXP_Ratio=[  8.7 9.3 10.8 12.3 16.8 20 ];
set(figure(4),'Units','inches','Position',[4 2 15 4]) 
figure(4)
subplot(1,3,1)
f4h1=plot(log10(dPO2),T_LAC,'k','LineWidth',line_width,'markersize',marker_size, 'MarkerFaceColor','k')

hold on 
f4h2= errorbar(log10(FEXP_PO2),FEXP_LAC,FEXP_LAC_SE,'d','marker','square','markersize',marker_size,...
                 'markeredgecolor','k', 'MarkerFaceColor','k',...
                'color','k','linewidth',2)
xlabel('log(PO_2 (mmHg))')
ylabel({'LAC production rate';'(\mumol/min/g dry weight)'})
legend([f4h1 f4h2],'Model','Data')
box off
% xlim([-0.1 50])
legend boxoff
%ylabel('Fractional ATP, ADP or AMP')
%title('LAC production rate')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
subplot(1,3,2)
plot(log10(dPO2),T_PYR,'k','LineWidth',line_width,'markersize',marker_size, 'MarkerFaceColor','k')
hold on

f4h3= errorbar(log10(FEXP_PO2),FEXP_PYR,FEXP_PYR_SE,'d','marker','square','markersize',marker_size,...
                 'markeredgecolor','k', 'MarkerFaceColor','k',...
                'color','k','linewidth',2)
            hold off
xlabel('log(PO_2(mmHg))')
ylabel({'PYR production rate';'(\mumol/min/g dry weight)'})
ylim([0 0.3])
% xlim([-0.1 50])
box off
legend boxoff
%ylabel('Fractional ATP, ADP or AMP')
%title('PYR production rate')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
subplot(1,3,3)
plot(log10(dPO2),T_LAC./T_PYR,'k','LineWidth',line_width,'markersize',marker_size, 'MarkerFaceColor','k')
xlabel('log(PO_2(mmHg))')
ylabel('LAC/PYR ratio')
hold on
h2s3= errorbar(log10(FEXP_PO2),FEXP_Ratio,FEXP_RATIO_SE,'d','marker','square','markersize',marker_size,...
                 'markeredgecolor','k', 'MarkerFaceColor','k',...
                'color','k','linewidth',2)
hold off
  box off
legend boxoff
%ylabel('Fractional ATP, ADP or AMP')
%title('LAC/PYR ratio')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)


set(figure(7),'Units','inches','Position',[4 2 15 4]) 
figure(7)
subplot(1,3,1)
plot(dPO2,T_LAC,'k','LineWidth',line_width)
hold on 
h2s1= errorbar(FEXP_PO2,FEXP_LAC,FEXP_LAC_SE,'d','marker','square','markersize',marker_size,...
                 'markeredgecolor','k', 'MarkerFaceColor','k',...
                'color','k','linewidth',2)
hold off
xlabel('PO_2 (mmHg)')
ylabel({'LAC production rate';'(\mumol/min/g dry weight)'})
legend('Model','Data')
box off
xlim([-0.1 800])
legend boxoff
%ylabel('Fractional ATP, ADP or AMP')
% title('LAC production rate')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)

subplot(1,3,2)
plot(dPO2,T_PYR,'k','LineWidth',line_width)
hold on 
h2s2= errorbar(FEXP_PO2,FEXP_PYR,FEXP_PYR_SE,'d','marker','square','markersize',marker_size,...
                 'markeredgecolor','k', 'MarkerFaceColor','k',...
                'color','k','linewidth',2)
            hold off
xlabel('PO_2 (mmHg)')
ylabel({'PYR production rate';'(\mumol/min/g dry weight)'})
xlim([-0.1 800])
ylim([0 0.3])
box off
%title('PYR production rate')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
subplot(1,3,3)
plot(dPO2,T_LAC./T_PYR,'k','LineWidth',line_width)
hold on
h2s3= errorbar(FEXP_PO2,FEXP_Ratio,FEXP_RATIO_SE,'d','marker','square','markersize',marker_size,...
                 'markeredgecolor','k', 'MarkerFaceColor','k',...
                'color','k','linewidth',2)
hold off
xlabel('PO_2 (mmHg)')
ylabel('LAC/PYR ratio')
xlim([-0.1 800])
ylim([0 30])
box off
%title('LAC/PYR ratio')
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)