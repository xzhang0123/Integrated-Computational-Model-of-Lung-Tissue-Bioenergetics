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

exp_glu05=[2.36 5.35 8.55 12.6 15.4]*Lung_dw;
exp_glu56=[10.8 23.5 35.7 56.2 69.5]*Lung_dw;
exp_glu15=[4.67 11 19.5 27.7 36.4]*Lung_dw;
exp_glu3=[7.2 16.3 28.8 41.9 53.5]*Lung_dw;
exp_glu10=[13.5 24.5 46.2 61.4 83.2]*Lung_dw;
%---------------dGLC------------
exp_lac05_exo=[13.1 24.8 30.4 32.6 33]*Lung_dw*0.18;  %Subtract endogenase part
exp_lac15_exo=[20.5 32.6 44.7 52.9 61]*Lung_dw*0.36;
exp_lac3_exo=[20.6 33.5 45.7 57.2 68.1]*Lung_dw*0.41;
exp_lac56_exo=[20 38.3 54.2 70.7 84.8]*Lung_dw*0.46;
exp_lac10_exo=[22.6 48.8 69 88.7 111]*Lung_dw*0.75;

exp_lac05=[13.1 24.8 30.4 32.6 33]*Lung_dw;  %Total
exp_lac15=[20.5 32.6 44.7 52.9 61]*Lung_dw; % ?????????????? lower than
exp_lac3=[20.6 33.5 45.7 57.2 68.1]*Lung_dw;
exp_lac56=[20 38.3 54.2 70.7 84.8]*Lung_dw;
exp_lac10=[22.6 48.8 69 88.7 111]*Lung_dw;

exp_lac05_endo=exp_lac05-exp_lac05_exo;  %Subtract endogenase part
exp_lac15_endo=exp_lac15-exp_lac15_exo; %
exp_lac3_endo=exp_lac3-exp_lac3_exo;
exp_lac56_endo=exp_lac56-exp_lac56_exo;
exp_lac10_endo=exp_lac10-exp_lac10_exo;
%----------------dLAC-------------------
exp_dLAC0=[16 32.7 53.4 75 95.4]*Lung_dw;
exp_dLAC05=[16 32.7 53.4 75 95.4]*Lung_dw;
exp_dLAC2=[7.41 15.8 26.9 37 49.3]*Lung_dw;
exp_dLAC5=[4.77 12.5 19.9 27.7 35.8]*Lung_dw;
exp_dLAC20=[3.11 9.54 16.6 24.4 31.5]*Lung_dw;
%------------------dPO2------------------------
    FEXP_LAC2=[71 69.8 86.3 85.1 112 132 ]/60*0.227; %unit umol/hr/g dry wt
      FEXP_PYR2=[7.9 7.4 7.9 6.7 6.6 6.6 ]/60*0.227;
  FEXP_Ratio2=[ 9.75 9.49 11.4 12.7 18.2 20.1];
%  dC_O2_step=0.003*FEXP_PO2/(22.4*1e6); %change Unit: M
 FEXP_LAC2=FEXP_LAC2/60*0.227;  %change unit to umol/min/lung
 %-------------------------------------------------------
size_index=dlmread('text_size.txt');
%-------------------ydata simulation2-----
% end of ydata simulation 2 ------------------------------------
 ydata=[exp_glu05 exp_glu15 exp_glu3 exp_glu56 exp_glu10 exp_lac05 exp_lac15 exp_lac3 exp_lac56 exp_lac10 ...
     exp_dLAC0 exp_dLAC05 exp_dLAC2 exp_dLAC5 exp_dLAC20 FEXP_LAC2 FEXP_PYR2  log10(FEXP_Ratio2)];
%
Para=ones(1,51);

Para_index=[1:9 13 49 46 32 33 34 45 48 ];
SSE0=Calculate_SSE(Para,ydata);
 for i=1:1:length(Para_index)
     
    Para1_loop=Para; Para2_loop=Para;
    Para1_loop(Para_index(i))=1.01*Para(Para_index(i));
    Para2_loop(Para_index(i))=0.99*Para(Para_index(i));

SSE1(i)=Calculate_SSE(Para1_loop,ydata);
SSE2(i)=Calculate_SSE(Para2_loop,ydata);
end
LSC=abs((SSE1-SSE2)./(0.02*SSE0));
LSC=LSC'