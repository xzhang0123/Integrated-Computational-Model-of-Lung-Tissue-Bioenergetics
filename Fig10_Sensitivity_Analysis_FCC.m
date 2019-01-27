clear all
clc
 
global Tem Flow F_con R_con Vr Vb Vc Vm Vi ...
   iGLCr iPYRr iLACr iPir iHr...
   iGLCb iPYRb iLACb iPib iHb...
   iGLCc iG6Pc iF6Pc iF16BPc iGAPc iBPGc iPEPc iPYRc iLACc iPG6c...
   iR5Pc iMALc iOXAc iCITc iaKGc iSUCc iFUMc iGLUc iASPc iPic...
   iAMPc iADPc iATPc iNADHc iNADc iNADPHc iNADPc iGSSG iGSH iH2O2...
    iHc   iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iGLUm...
    iASPm iNADm iNADHm iACOAm iCOAm iUQm iUQH2m  iPim iADPm iATPm...
    iFADm iFADH2m iHm iCytCoxi iCytCred iHi idPsim idPsip iO2 iR123e iR123m

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
Lung_dw=0.227;
  %------------
    Flow   =   12e-3;     %mL/min  
    %LUNG WEIGH 1.6 mL,  Vb is 10%  Vc is 90%, Vmito=~2% cells Vi=~10% of Mit
Vr  =   55e-3;      %mL
%Lung_volume=1.33e-3;    %1.6mL 
Vb  =   0.66e-3; %mL
Vcell  =  0.67e-3; %mL  

Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]

  Tem=310.15; %K      37 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]

% Inhibitor_factor=[ 1 :-0.1:0.1    0.06   0.05];  %
Inhibitor_factor=[ 1 :-0.05:0.2]; 
 %-------------------------------------------------------
size_index=dlmread('text_size.txt');
%-------------------ydata simulation2-----
% end of ydata simulation 2 ------------------------------------


for loop=1:1:length(Inhibitor_factor)
Para=ones(1,51);
Para(27)=Inhibitor_factor(loop)*Para(27);%CI inhibitor
Para_index=[1 :8 32];
%Para_index=[1 3 7];
[EC,SSE0]=Calculate_FCC(Para);
 for i=1:1:length(Para_index)
     
    Para1_loop=Para; Para2_loop=Para;
    Para1_loop(Para_index(i))=1.01*Para(Para_index(i));
    Para2_loop(Para_index(i))=0.99*Para(Para_index(i));

[~,SSE1(i)]=Calculate_FCC(Para1_loop);
[~,SSE2(i)]=Calculate_FCC(Para2_loop);
end
LSC(:,loop)=abs((SSE1-SSE2)./(0.02*SSE0));
ALL_Other(loop)=sum(LSC([2 4:6 8],loop));
%LSC=LSC'
ECi(loop)=EC;
 end

 figure(1)
 plot(Inhibitor_factor,LSC)
 xlabel('CI activity')
ylabel('Flux control coefficent')
title('Flux control coefficent')
box off
 legend('HK','PGI','PFK','ALD','GAPDH','PGK','PK','LDH','GLUT')
 legend boxoff
   set(gcf,'color','w')
set(gca,'Fontsize',16)

 figure(2)
 plot(ECi,LSC)
 xlabel('Energy charge')
ylabel('Flux control coefficent')
box off
 legend('HK','PGI','PFK','ALD','GAPDH','PGK','PK','LDH','GLUT')
 legend boxoff
   set(gcf,'color','w')
set(gca,'Fontsize',16)

set(figure(3),'Units','inches','Position',[0.2 0.1 5 6]) 
set(figure(4),'Units','inches','Position',[0.2 0.1 5 4]) 
 figure(3)
 plot(ECi,LSC(1,:),'r',ECi,LSC(3,:),'g',ECi,LSC(7,:),'-b',ECi,LSC(9,:),'--k','LineWidth',2)
 xlabel('Energy charge')
ylabel('Flux control coefficent')
%title('Flux control coefficent')
box off
 legend('HK','PFK','PK','GLUT')
 legend boxoff
   set(gcf,'color','w')
set(gca,'Fontsize',16)



 figure(4)
 plot(Inhibitor_factor,LSC(1,:),'r',Inhibitor_factor,LSC(3,:),'g',Inhibitor_factor,LSC(7,:),'b',Inhibitor_factor,LSC(9,:),'m','LineWidth',2)
 xlabel('%CI')
ylabel('Flux control coefficent')
title('Flux control coefficent')
box off
 legend('HK','PFK','PK','GLUT')
 legend boxoff
   set(gcf,'color','w')
set(gca,'Fontsize',16)
% 
%  figure(5)
%  plot(ECi,sum(LSC))
%  xlabel('Energy charge')
% ylabel('Flux control coefficent')
% title('Flux control coefficent')
% box off
%  legend('HK','PFK','PK','GLUT','Sum of other reactions')
%  legend boxoff
%    set(gcf,'color','w')
% set(gca,'Fontsize',16)