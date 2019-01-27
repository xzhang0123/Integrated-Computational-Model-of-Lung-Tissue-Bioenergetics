function rtfluxes=fluxes(C_t,Para)


 global F_con R_con Tem Vr Vb Vc Vm Vi ...

   
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
    
    r1=1;r2=2;r3=3;r4=4;r5=5;r6=6;r7=7;r8=8;r9=9;r10=10;r11=11;r12=12;r13=13;
    r14=14;r15=15;r16=16;r17=17;r18=18;r19=19;r20=20;r21=21;r22=22;r23=23;r24=24;
    r25=25;r26=26;r27=27;r28=28;r29=29;r30=30;r31=31;r32=32;
    
    Mass_mito=16.443564787283922*Para(48); %Mass of mitochondria

    %----initialization
    Vmaxf=1e-6*ones(1,31);
    Tmax=1e-6*ones(1,16);

   %----Set parameters:glycolysis----------------------
    Vmaxf(1)=2.790450000000000e+02*1e-9;      %HEX
    Vmaxf(2)=3.701880000000000e+04*1e-9;      %PGI
    Vmaxf(3)=0.225*1.467214320156125e+03*1e-9;      %PFK
    Vmaxf(4)=1.426207679135680e+03*1e-9;      %ALD
    Vmaxf(5)=6.349528210139666e+04*1e-9;      %G3PDH
    Vmaxf(6)=2.002726365585151e+05*1e-9;      %PHK
    Vmaxf(7)=1.888287076523587e+04*1e-9;      %PK
    Vmaxf(8)=3.919607740112311e+05*1e-9;      %LDH
    %----PPP pathway-------------------------
    Vmaxf(9)=2.826212399999999e+03*1e-9;  %PPP1  Para not known 
    Vmaxf(10)=0;    %PPP2   not used
    %----------------------
    Vmaxf(11)=7.762499999999999e-4;                %GSH1
    Vmaxf(12)=9e-6;                %GSH2
    Vmaxf(13)=2.536261317373874e-05;          %ATPase 
    Vmaxf(14)=10e-3;                          %AK  This reaction is at equalibrium
    Vmaxf(15:16)=0;  % not used, reserved for future use
    %------Mitochondria reaction parameters-------------------
Vmaxf(17:31)=[
    306.49e-09;                               %PDH
    4630.8*1e-9;                              %CITS
    530.797548526789e-09;                     %CITDH
    780e-9;                                   %AKGDH 
    5.81695925037748e5*1e-9;                  %SCAS
    4.32513873910319e5*1e-09;                 %NDK  equalibrium reaction
    1.08593076185781e-05;                     %SDH
    6435e-9;                                  %FH  
    3297.559540729090*1e-9;                   %MDH
   	975e-9;                                   %got,value obtained from ref=Sources of C for TCA
    11.007e-9;                                %CI
    220*1e-9;                                 %CII  ref=Age related changes in metabolism oxidative stress
    2.2335e4*1e-9;                            %CIII
    0.271043306907577*1e-9;                   %CIV
    588.98e-09; ]  ;                          %CV
%-------------------------------------
%Unit: nmol/min/lung at 30 degree
    Tmax(1)=5.341245681509842e+02*1e-9;       %GLC 
    Tmax(2)=1.850877922494081e+02*1e-9;       %PYR
    Tmax(3)=0.7*0.9*1.009745540579863e+04*1e-9;       %LAC
    Tmax(4)=1.2e-3*100*1e-6;                         %PI (bc), equlibrium, set to a big value 
    %---------------------------
    Tmax(5)=1698.5e-09;                       %SUC PI
    Tmax(6)= 21.98e-09;                       %MAL PI
    Tmax(7)= 0.9994e-09;                      %OME 
    Tmax(8)=81.33e-09;                        %TCC 
    Tmax(9)= 96.62e-09;                       %PYRH 
    Tmax(10)=1*2.38e-5;                         %PIC
    Tmax(11)=523.8837274961545e-09;           %ANT
    Tmax(12)=5*1e-9;                          %GLUH
    Tmax(13)=646*1e-9;                        %ASP GLU
    Tmax(14)=0.45*35.5*1e-9;                  %LEAK Para 45
    Tmax(15)=4.428614065699552e+03*1e-9;      %Lumped MA shuttle
    Tmax(16)=0;                               %ALA reaction
    
   Tmax(5:7)=0;  %LUMPED to flux T15
   Tmax(12:13)=0;    %Lumped to MAS
 %_--Adjust parameters (used for parameter estimation)------
 Vmaxf=Vmaxf.*Para(1:31);
 Tmax=Tmax.*Para(32:47);
%--- Change unit from (Mole/mass of mito) to(Mole/lung)
   Vmaxf(17:31)=Mass_mito.*Vmaxf(17:31);
   Tmax(5:14)=Mass_mito.*Tmax(5:14);
 %--------------Temperature correction--------
Q10=2.5;
Tem_standard=303.15;  %Default temperature (30oC), current Tem=37
Q10_factor=Q10.^((Tem-Tem_standard)/10);
Tmax=Tmax.*Q10_factor;
Vmaxf=Vmaxf.*Q10_factor;
%-----------------------------------------
RT=R_con*Tem;
%%  Thermodymamic data (Gibbs free energy of reactions, deltaG)
dGr(1)= -24.41;  %HEX
dGr(2)= 3.18;    %PISO
dGr(3)= -23.24;    %PFK
dGr(4)= 30.7;      %G3PF
dGr(5)= 1.21;     %G3PD
dGr(6)= -5.85;     %PHK
dGr(7)= -28.89;   %PK
dGr(8)= -23.9;   %LD
dGr(9)= -54.89;   %PPP1 
dGr(10)= 1.23; %PPP2   not used
dGr(11)= 0;   %GSR Keq=0.075, defined later
dGr(12)= -275.7;   %GP
dGr(13)= - 36.03;  %ATPase
dGr(14)= -0.475;    %AK
dGr(15)= -28.83;     %MD
dGr(16)= 1.4;     %GOT2
dGr(17)= -38.64; %PYR dehydrogenase [1]
dGr(18)= -36.61;%CIT synthase [1]
dGr(19)=2.18;%aconiase+isocit dehydrogenase [1]
dGr(20)=-37.08;%aKG dehydro
dGr(21)=1.26;%scoa synthetase
dGr(22)=0;   % GDP+ATP-GTP-ADP
dGr(23)=-3.62;%suc dehydro
dGr(24)=-3.6; %FUM Jason's model
dGr(25)=28.83;% MDH2
dGr(26)=-1.4;%asp+akg-glu+oxa [3]
dGr(27)=-69.37;%CI   [2] 
dGr(28)= -1.31;   %[2]-[1]
dGr(29)=-32.53; %CIII  [2]
dGr(30)=-122.94;%CIV   [2]
dGr(31)=36.03;%CV [2]
dGr(32)=-20;
%referneces:
%[1] Role of NADH/NAD transport 
%[2] analysis of cardiac mito NA/Ca exchange kinetics with a biophysical
%model of ~~~~ 3:1 stoichiometry
%[3] Appendices for ?Computer Modeling of Mitochondrial TCA Cycle, Oxidative
%Phosphorylation, Metabolite Transport, and Electrophysiology?
Keq0=exp(-dGr/RT);
Keq0(11)=0.075;  %GSR
CO2=1.32e-5;
%-----------
GLCr=C_t(iGLCr);
PYRr=C_t(iPYRr);
LACr=C_t(iLACr);
Pir=C_t(iPir);
Hr=C_t(iHr);
%--------------------
GLCb=C_t(iGLCb);
PYRb=C_t(iPYRb);
LACb=C_t(iLACb);
Pib=C_t(iPib);
Hb=C_t(iHb);
%----------------------
GLCc=C_t(iGLCc);
G6Pc=C_t(iG6Pc);
F6Pc=C_t(iF6Pc);
F16BPc=C_t(iF16BPc);
GAPc=C_t(iGAPc);
BPGc=C_t(iBPGc);
PEPc=C_t(iPEPc);
PYRc=C_t(iPYRc);
LACc=C_t(iLACc);
PG6c=C_t(iPG6c);
R5Pc=C_t(iR5Pc);
MALc=C_t(iMALc);
OXAc=C_t(iOXAc);
CITc=C_t(iCITc);
aKGc=C_t(iaKGc);
SUCc=C_t(iSUCc);
FUMc=C_t(iFUMc);
GLUc=C_t(iGLUc);
ASPc=C_t(iASPc);
Pic=C_t(iPic);
AMPc=C_t(iAMPc);
ADPc=C_t(iADPc);
ATPc=C_t(iATPc);
NADHc=C_t(iNADHc);
NADc=C_t(iNADc);
NADPHc=C_t(iNADPHc);
NADPc=C_t(iNADPc);
GSSGc=C_t(iGSSG);
GSHc=C_t(iGSH);
H2O2=C_t(iH2O2);
Hc=C_t(iHc);
%------------------
PYRm = C_t(iPYRm);
OXAm = C_t(iOXAm);
CITm = C_t(iCITm);
aKGm = C_t(iaKGm);
SCAm = C_t(iSCAm);
SUCm = C_t(iSUCm);
FUMm = C_t(iFUMm);
MALm = C_t(iMALm);
GLUm = C_t(iGLUm);
ASPm = C_t(iASPm);
NADm  = C_t(iNADm);
NADHm = C_t(iNADHm);
ACOAm= C_t(iACOAm);
COAm= C_t(iCOAm);
UQm   = C_t(iUQm);
UQH2m =C_t(iUQH2m);
Pim = C_t(iPim);
ADPm = C_t(iADPm);
ATPm = C_t(iATPm);
GDPm=ADPm;
GTPm=ATPm;
FADm =C_t(iFADm);
FADH2m = C_t(iFADH2m);    
Hm=C_t(iHm);
%------------------
CytCoxi=C_t(iCytCoxi);
CytCred=C_t(iCytCred);
Hi=C_t(iHi);
dPsim=C_t(idPsim);
dPsip=C_t(idPsip);
O2 =C_t(iO2);
R123e=C_t(iR123e);
R123m=C_t(iR123m);
DeltaGH = real(F_con*dPsim+ R_con*Tem*log(Hi/Hm));
pH_c=-log10(Hc);
pH_m=-log10(Hm);

%----------------------------------------------
%% Chemical Reaction Rate Equations, R(n), by Vmax(n)
%R1: GLCc+ATP->G6P+ADP+Hc
KA(r1)=1e-4;
KB(r1)=1.16e-4;
KC(r1)=1e-4;
KD(r1)=1.26e-4;
A=GLCc;
B=ATPc;
C=G6Pc;
D=ADPc;
Keq(r1)=Keq0(r1)*10^(pH_c-7);
deno=(1+A/KA(r1))*(1+B/KB(r1)+D/KD(r1))*(1+C/KC(r1));
Ki_HEX=27.8e-3;   %value obtained from ref : NADH/NAD shuttle
%-------------------------
Ctrl_HEX=Ki_HEX/(Ki_HEX+G6Pc);
J_HEX=Ctrl_HEX*Vmaxf(r1)/KA(r1)/KB(r1)*(A*B-C*D/Keq(r1))/deno;
%----------------------------------------------
%R2: G6P->F6P
KA(r2)=0.59e-3;
KB(r2)=0.095e-3;
A=G6Pc;
B=F6Pc;
Keq(r2)=Keq0(r2);
deno=(1+A/KA(r2))*(1+B/KB(r2));
J_PGI=Vmaxf(r2)/KA(r2)*(A-B/Keq(r2))/deno;
%----------------------------------------------
%R3: F6P+ATP->F16BP+ ADP+H+  (pfk)
KA(r3)=32e-6;
KB(r3)=20e-6;
KC(r3)=32e-6;
KD(r3)=20e-6;
A=F6Pc;
B=ATPc;
C=F16BPc;
D=ADPc;
Keq(r3)=Keq0(r3)*10^(pH_c-7);
deno=(1+A/KA(r3))*(1+B/KB(r3)+D/KD(r3))*(1+C/KC(r3));
% % %-------------------Cooperative control mechansim----------
hcef=2.3;
Ki_ATP=0.072*0.068e-3; %
Ki_AMP=1.2*0.05*0.035e-3;%
Ki_ADP=0.9*0.12*0.035e-3;%
Ki_CIT=2.6e-3;
L_factor=(1+CITc/Ki_CIT).^(hcef)*(1+ATPc/Ki_ATP).^(hcef)./(1+F6Pc/KA(r3)+F16BPc/KC(r3)).^hcef./(1+AMPc/Ki_AMP)^(hcef)./(1+ADPc/Ki_ADP)^(hcef);
Ctrl_PFK=1./(1+L_factor);
J_PFK=Ctrl_PFK*Vmaxf(r3)/KA(r3)/KB(r3)*(A*B-C*D/Keq(r3))/deno;
% %----------------------------------------------
%R4: F16BP-> 2GA3P   (ALD)
KA(r4)=1.7e-6;
KB(r4)=1.7e-6;
A=F16BPc;
B=GAPc;
Keq(r4)=Keq0(r4);
deno=(1+A/KA(r4))*(1+B^2/(KB(r4).^2));
J_G3PF=Vmaxf(r4)/KA(r4)*(A-B.^2/Keq(r4))/deno;
%----------------------------------------------
%R5: GAP+Pi+NADc->13BPG+NADH +Hc
KA(r5)=20.7e-6;
KB(r5)=2.9e-3;  %set to intitial concentration
KC(r5)=17.8e-6;
KD(r5)=20.7e-6;
KE(r5)=17.8e-6;
A=GAPc;
B=Pic;
C=NADc;
D=BPGc;
E=NADHc;
Keq(r5)=Keq0(r5)*10^(pH_c-7);
deno=(1+A/KA(r5))*(1+B/KB(r5))*(1+C/KC(r5)+E/KE(r5))*(1+D/KD(r5));
J_G3PD=Vmaxf(r5)/KA(r5)/KB(r5)/KC(r5)*(A*B*C-D*E/Keq(r5))/deno;
%----------------------------------------------
%R6: BPG+ADP->PEP+ATP   (phk)
KA(r6)=10e-6;
KB(r6)=0.4e-3;
KC(r6)=10e-6;
KD(r6)=0.4e-3;
A=BPGc;
B=ADPc;
C=PEPc;
D=ATPc;
Keq(r6)=Keq0(r6);
deno=(1+A/KA(r6))*(1+B/KB(r6)+D/KD(r6))*(1+C/KC(r6));
J_PHK=Vmaxf(r6)/KA(r6)/KB(r6)*(A*B-C*D/Keq(r6))/deno;
%----------------------------------------------
%R7: PEP+ADP+Hc->PYR+ATP  (pk)
KA(r7)=0.26e-3;
KB(r7)=0.4e-3;
KC(r7)=0.26e-3;
KD(r7)=0.4e-3;
A=PEPc;
B=ADPc;
C=PYRc;
D=ATPc;
Keq(r7)=Keq0(r7)*10^(7-pH_c);

%------------------------------------
Ki_PK_F16BP=1e-6;  %value obtained from ref
Ctrl_PK=F16BPc/(Ki_PK_F16BP+F16BPc);
%------------------------------------
% Ctrl_PK=0.5;
deno=(1+A/KA(r7))*(1+B/KB(r7)+D/KD(r7))*(1+C/KC(r7));
J_PK=Ctrl_PK*Vmaxf(r7)/KA(r7)/KB(r7)*(A*B-C*D/Keq(r7))/deno;
%----------------------------------------------
%R8: PYR+NADH+Hc->LAC+NAD  (LD)
KA(r8)=0.052e-3;  %PYR
KB(r8)=0.68e-6;   %NADH
KC(r8)=4.934e-3;
KD(r8)=340e-6;
A=PYRc;
B=NADHc;
C=LACc;
D=NADc;
Keq(r8)=Keq0(r8)*10^(7-pH_c);
deno=(1+A/KA(r8)+C/KC(r8))*(1+B/KB(r8)+D/KD(r8));
J_LD=Vmaxf(r8)/KA(r8)/KB(r8)*(A*B-C*D/Keq(r8))/deno;
%----------------------------------------------
 %R9&10  1/6 G6P+2 NADP->CO2+2 NADPH + 2Hc     lumped PPP reaction
KA(r9)=36e-6;
KB(r9)=4.8e-6;
KC(r9)=36e-6;
KD(r9)=4.8e-6;
A=G6Pc;
B=NADPc;
C=CO2;
D=NADPHc;
Keq(r9)=Keq0(r9)*10^(pH_c-7)*10^(pH_c-7);
deno=(1+(A/KA(r9))^(1/6))*(1+(B/KB(r9))^2+(D/KD(r9))^2)*(1+C/KC(r9));
J_PPP1=Vmaxf(r9)/KA(r9).^(1/6)/KB(r9).^2*(A.^(1/6)*B.^2-C*D.^2/Keq(r9))/deno;
J_PPP2=0; %not used 
%----------------------------------------------
%R11: GSSG+NADPH+Hc->2GSH+NADP (GSR)
KA(r11)=61e-6;
KB(r11)=7.6e-6;
KC(r11)=61e-6;
KD(r11)=7.6e-6;
A=GSSGc;
B=NADPHc;
C=GSHc;
D=NADPc;
Keq(r11)=Keq0(r11)*10^(7-pH_c);
deno=(1+A/KA(r11))*(1+B/KB(r11)+D/KD(r11))*(1+(C/KC(r11))^2)*(1);
J_GSR=Vmaxf(r11)/KA(r11)/KB(r11)*(A*B-C*C*D/Keq(r11))/deno;
%----------------------------------------------
%R12: 2GSH+H2O2->GSSG+2H2O (GP)
KA(r12)=3e-3;
KB(r12)=1e-6;
KC(r12)=3e-3;
A=GSHc;
B=H2O2;
C=GSSGc;
Keq(r12)=Keq0(r12); 
deno=(1+(A/KA(r12))^2+C/KC(r12))*(1+B/KB(r12));
J_GP=Vmaxf(r12)/KA(r12)/KB(r12)*(A*A*B-C/Keq(r12))/deno;
%----------------------------------------------
%R13: ATP->ADP+Pi+HC
KA(r13)=1e-3;
KB(r13)=100e-6;
KC(r13)=1e-3;
A=ATPc;
B=ADPc;
C=Pic;
Keq(r13)=Keq0(r13)*10^(pH_c-7);
deno=(1+A/KA(r13))*(1+B/KB(r13))*(1+C/KC(r13));
J_ATPASE=Vmaxf(r13)/KA(r13)*(A-B*C/Keq(r13))/deno;
%----------------------------------------------
%R14: AMP+ATP->2ADP   AK
KA(r14)=10e-6;
KB(r14)=1e-3;
KC(r14)=100e-6;
A=AMPc;
B=ATPc;
C=ADPc;

Keq(r14)=Keq0(r14);
deno=(1+A/KA(r14))*(1+B/KB(r14))*(1+C/KC(r14))*(1+C/KC(r14));
J_AK=Vmaxf(r14)/KA(r14)/KB(r14)*(A*B-C*C/Keq(r14))/deno;
%----------------------------------------------
%R15: OXAc+NADH+HC--->MAL+NAD   reverse MDH
KA(r15)=3.38e-6;
KB(r15)=3.61e-5;
KC(r15)=1.55e-4;
KD(r15)=1.1e-3;
A=OXAc;
B=NADHc;
C=MALc;
D=NADc;
Keq(r15)=Keq0(r15)*10^(7-pH_c);
deno=(1+A/KA(r15))*(1+B/KB(r15))*(1+C/KC(r15))*(1+D/KD(r15));
J_MD2=Vmaxf(r15)/KA(r15)/KB(r15)*(A*B-C*D/Keq(r15))/deno;
%----------------------------------------------
%R16: GLUc+OXA->ASP+aKGc     reversed GOT
KA(r16)=8.9e-3;
KB(r16)=88e-6;
KC(r16)=3.9e-3;
KD(r16)=430e-6;
A=GLUc;
B=OXAc;
C=ASPc;
D=aKGc;
Keq(r16)=Keq0(r16);
deno=(1+A/KA(r16))*(1+B/KB(r16))*(1+C/KC(r16))*(1+D/KD(r16));

J_GOT2=Vmaxf(r16)/KA(r16)/KB(r16)*(A*B-C*D/Keq(r16))/deno;
%----------------------------------------------
%R17: PYR+COAm+NAD->ACOA+CO2+NADH+Hm
KA(r17)=5.08e-3;
KB(r17)=1.49e-5;
KC(r17)=3.5e-5;
KD(r17)=1.49e-5;
KE(r17)=1e-5;%not used in model
KF(r17)=3.5e-5;
A=PYRm;
B=COAm;
C=NADm;
D=ACOAm;
E=CO2;
F=NADHm;
Keq(r17)=Keq0(r17)*10^(pH_m-7);
deno=(1+A/KA(r17))*(1+B/KB(r17)+D/KD(r17))*(1+C/KC(r17)+F/KF(r17));
J_PDH = Vmaxf(r17)/KA(r17)/KB(r17)/KC(r17)*(A*B*C-D*E*F/Keq(r17))/deno;
%-----------------------------------------------
%R18: ACOA+OXA->CIT+COA+Hm
KA(r18)=3.9e-6;
KB(r18)=4.53e-6;
KC(r18)=57.9e-6;
KD(r18)=4.3e-3;
A=ACOAm;
B=OXAm;
C=COAm;
D=CITm;
Keq(r18)=Keq0(r18)*10^(pH_m-7)*10^(pH_m-7);
deno=(1+A/KA(r18)+C/KC(r18))*(1+B/KB(r18))*(1+D/KD(r18));
J_CITS = Vmaxf(r18)/KA(r18)/KB(r18)*(A*B-C*D/Keq(r18))/deno;
%----------------------------------------------
%R19: CIT+NAD+->aKGm+NADHm+CO3--+Hm
KA(r19)=1.3e-3;
KB(r19)=500e-6;
KC(r19)=3.5e-6;
KD(r19)=4.7e-6;
KE(r19)=1e-5; %CO2,NOT USED
A=CITm;
B=NADm;
C=aKGm;
D=NADHm;
E=CO2;
Keq(r19)=Keq0(r19);
deno=(1+A/KA(r19))*(1+C/KC(r19))*(1+B/KB(r19)+D/KD(r19));
J_CITD = Vmaxf(r19)/KA(r19)/KB(r19)*(A*B-C*D*E/Keq(r19))/deno;
%----------------------------------------------
%R20: aKG + CoA+NAD+H2O-> SCoA+ NADH+CO3+H+
KA(r20)=85.87e-6;
KB(r20)=1.634e-6;
KC(r20)=46.6e-6;
KD(r20)=KB(r20);
KE(r20)=7.14e-6;
KF(r20)=1e-5;%CO2,NOT USED
A=aKGm;
B=COAm;
C=NADm;
D=SCAm;
E=NADHm;
F=CO2;
Keq(r20)=Keq0(r20)*10^(pH_m-7);
deno=(1+A/KA(r20)+D/KD(r20))*(1+B/KB(r20)+E/KE(r20))*(1+C/KC(r20));
J_AKGD =Vmaxf(r20)/KA(r20)/KB(r20)/KC(r20)*(A*B*C-D*E*F/Keq(r20))/deno;
%----------------------------------------------
%R21 SCA+GDPm+Pim ->SUC+GTP+COA+H+
KA(r21)=1.2538e-5;
KB(r21)=5.6977e-6;
KC(r21)=2.1e-3;
KD(r21)=4.51e-4;
KE(r21)=2.3143e-5;
KF(r21)=1.7735e-5;
A=SCAm;
B=GDPm;
C=Pim;
D=SUCm;
E=GTPm;
F=COAm;
Keq(r21)=Keq0(r21)*10^(pH_m-7);
deno=(1+A/KA(r21)+D/KD(r21)+E/KE(r21)+D*E/KD(r21)/KE(r21))*(1+B/KB(r21)+C/KC(r21)+F/KF(r21)+B*C/KB(r21)/KC(r21));
J_SCAS= Vmaxf(r21)/KA(r21)/KB(r21)/KC(r21)*(A*B*C-D*E*F/Keq(r21))/deno;
%----------------------------------------------
%R22 GTP+ADP ----GDP+ATP
KA(r22)=111e-6;
KB(r22)=100e-6;
KC(r22)=260e-6;
KD(r22)=278e-6;
A=GTPm;
B=ADPm;
C=GDPm;
D=ATPm;
Keq(r22)=Keq0(r22);
deno=(1+A/KA(r22)+C/KC(r22))*(1+B/KB(r22)+D/KD(r22));
J_NDK=Vmaxf(r22)/KA(r22)/KB(r22)*(A*B-C*D/Keq(r22))/deno;
%----------------------------------------------
%R23 SUC+FAD->FUM+FADH2
KA(r23)=1800e-6;
KB(r23)=140e-6;
KC(r23)=1800e-6;
KD(r23)=2.45e-6;
A=SUCm;
B=FADm;
C=MALm;
D=FADH2m;
Keq(r23)=Keq0(r23);
deno=(1+A/KA(r23))*(1+C/KC(r23))*(1+B/KB(r23)+D/KD(r23));
J_SUCD =Vmaxf(r23)/KA(r23)/KB(r23)*(A*B-C*D/Keq(r23))/deno;
%-----------------------------------------------
%R24 FUM+H2O->MAL
KA(r24)=1800e-6;
KB(r24)=1800e-6;
A=FUMm;
B=MALm;
Keq(r24)=Keq0(r24);
deno=(1+A/KA(r24))*(1+B/KB(r24));
J_FUM =Vmaxf(r24)/KA(r24)/KB(r24)*(A-B/Keq(r24))/deno;
%----------------------------------------------
%R25 MAL+NAD+->OXA+NADH+Hm
KA(r25)=1.55e-4;
KB(r25)=1.1e-3;
KC(r25)=3.38e-6;
KD(r25)=3.61e-5;
A=MALm;
B=NADm;
C=OXAm;
D=NADHm;
Keq(r25)=Keq0(r25)*10^(pH_m-7);
deno=(1+A/KA(r25)+C/KC(r25))*(1+B/KB(r25)+D/KD(r25));
J_MALD =Vmaxf(r25)/KA(r25)/KB(r25)*(A*B-C*D/Keq(9))/deno;
%----------------------------------------------
% R26 ASPm+aKGm-GLUm+OXA
KA(r26)=3.9e-3;  %no data on this reaction
KB(r26)=430e-6;
KC(r26)=8.9e-3;
KD(r26)=88e-6;
A=ASPm;
B=aKGm;
C=GLUm;
D=OXAm;
Keq(r26)=Keq0(r26);
deno=(1+A/KA(r26))*(1+B/KB(r26))*(1+C/KC(r26))*(1+D/KD(r26));
J_GOT =Vmaxf(r26)/KA(r26)/KB(r26)*(A*B-C*D/Keq(r26))/deno;
%----------------------------------------------
%R11 Complex I   NADHm+UQm 5Hm+--NAD+UQH2+4Hi+
betaC1=0.5;
KA(r27)=1.5e-6;
KB(r27)=58.1e-6;
KC(r27)=428e-6;
KD(r27)=520e-6;
A=NADHm;
B=UQm;
C=NADm;
D=UQH2m;
Keq(r27)=Keq0(r27)*10^(7-pH_m);

deno=(1+A/KA(r27)+C/KC(r27))*(1+B/KB(r27)+D/KD(r27));
a=Keq(r27).^0.5;
J_CI =Vmaxf(r27)/KA(r27)/KB(r27)*(a*exp(-4*betaC1*DeltaGH/RT)*A*B-exp(-4*(betaC1-1)*DeltaGH/RT)*C*D/a)/deno;
%----------------------------------------------
% Complex II    FADH2+UQm--FADm+UQH2m
KA(r28)=1.5e-6;
KB(r28)=58e-6;
KC(r28)=1.5*285e-6;
KD(r28)=4*130e-6;
A=FADH2m;
B=UQm;
C=FADm;
D=UQH2m;
Keq(r28)=Keq0(r28);
deno=(1+A/KA(r28)+C/KC(r28))*(1+B/KB(r28)+D/KD(r28));
J_CII=Vmaxf(r28)/KA(r28)/KB(r28)*(A*B-C*D/Keq(r28))/deno;
%----------------------------------------------
% Complex III   UQH2m+2CytCoxi+2Hm--UQm+2CytCred+4Hi
betaC3=0.5;
KA(r29)=4.66e-6;
KB(r29)=3.76e-6;
KC(r29)=4.08e-6;
KD(r29)=4.91e-6;
A=UQH2m;
B=CytCoxi;
C=UQm;
D=CytCred;
Keq(r29)=Keq0(r29)*10^(pH_m-7)*10^(pH_m-7);%
deno=(1+A/KA(r29)+C/KC(r29))*(1+B^2/KB(r29)^2+D^2/KD(r29)^2);
J_CIII=Vmaxf(r29)/KA(r29)/(KB(r29)^2)*(Keq(r29).^(0.5)*exp(betaC3*(-4*DeltaGH/RT+2*F_con*dPsim/RT))*A*B^2-exp((betaC3-1)*(-4*DeltaGH/RT+2*F_con*dPsim/RT))*C*D^2*Keq(r29).^(-0.5))/deno;
%----------------------------------------------
%Complex IV     2CytCred+0.5 O2+4Hm---2CytCoxi+H2O+2Hi
betaC4=0.5;
KA(r30)=680e-6;
KB(r30)=5.4e-6*Para(51);
KC(r30)=79.2e-6;
A=CytCred;
B=O2;
C=CytCoxi;
Keq(r30)=Keq0(r30)*10^(7-pH_m)*10^(7-pH_m);
deno=(1+A^2/KA(r30)^2+C^2/KC(r30)^2)*(1+B^0.5/KB(r30)^0.5);
J_CIV =Vmaxf(r30)/KA(r30)^2/KB(r30)^0.5*(Keq(r30).^0.5*exp(betaC4*(-2*DeltaGH/RT-2*F_con*dPsim/RT))*A^2*B^0.5-exp((betaC4-1)*(-2*DeltaGH/RT-2*F_con*dPsim/RT))*C^2*Keq(r30).^(-0.5))/deno;
%----------------------------------------------
%ADPm+Pi+3Hi+Hm->ATPm+3Hm
betaC5=0.5;
KA(r31)=1e-3;% no data
KB(r31)=3e-3;
KC(r31)=1e-3;
A=ADPm;
B=Pim;
C=ATPm;
Keq(r31)=Keq0(r31)*10^(7-pH_m);
deno=(1+A/KA(r31)+C/KC(r31))*(1+B/KB(r31));
J_CV =Vmaxf(r31)/KA(r31)/KB(r31)*(Keq(r31).^0.5*exp(3*betaC5*DeltaGH/RT)*A*B-exp(3*(betaC5-1)*DeltaGH/RT)*C*Keq(r31).^(-0.5))/deno;
%--ALA--Transamination
 Ki_ALA_ANP=0.0874; %initial ratio
Ctrl_ALA=(AMPc/ATPc)/(Ki_ALA_ANP+AMPc/ATPc);  %
%--------------------------------------
K_PYR=1/40*1e-3;
K_ALA=1.3e-3; % ALA concentration constant due to degradation of  protein
ALAc=5e-3;
deno=(1+PYRc/K_PYR)*(1+ALAc/K_ALA);
Vmax_ALA=40.9e2*1e-09*Para(49);
Keq(r32)=Keq0(r32);
J_ALA=Ctrl_ALA*Vmax_ALA/K_ALA*(ALAc-PYRc/Keq(r32))/deno;

%% Transport fluxess

%between b and c 
%------------------------
 Ki_GLUT_ANP=0.077173815917969*Para(50); %initial ratio
Ctrl_GLUT=(AMPc/ATPc)/(Ki_GLUT_ANP+AMPc/ATPc);
T_GLCT_bc=Ctrl_GLUT*Tmax(1)*(GLCb-GLCc)/(3.4e-3+ (GLCb+GLCc));
%--------------------
 KPYRTbc_PYR=0.1e-3;      %ref= "Transpulmonary pyruvate kinetics"
KPYRTbc_H=1e-7;              
deno=1+PYRb/KPYRTbc_PYR+PYRc/KPYRTbc_PYR+Hb/KPYRTbc_H+Hc/KPYRTbc_H+PYRb*Hb/KPYRTbc_H/KPYRTbc_PYR+PYRc*Hc/KPYRTbc_H/KPYRTbc_PYR;
T_PYRT_bc= Tmax(2)/KPYRTbc_PYR/KPYRTbc_H*(Hb*PYRb-Hc*PYRc)/deno;
%--------------------
KLACTbc_LAC=4e-3;       %4mM, ref= "Transpulmonary pyruvate kinetics"
KLACTbc_H=1e-7;
deno=1+LACb/KLACTbc_LAC+LACc/KLACTbc_LAC+Hb/KLACTbc_H+Hc/KLACTbc_H+LACb*Hb/KLACTbc_H/KLACTbc_LAC+LACc*Hc/KLACTbc_H/KLACTbc_LAC;
T_LACT_bc=Tmax(3)/KLACTbc_LAC/KLACTbc_H*(LACb*Hb-LACc*Hc)/deno;
%--------------------
T_PIT_bc=Tmax(4)*(Pib-Pic)/(2.9e-3+Pib +Pic);
%-------------------------------------------------
%SUC--Pi   SUCe+Pim -- SUCm+Pie   Antipoter   A1+B2=A2+B1
KDCC_PI=0.93e-3;
KDCC_SUC=1e-3; 
KDCC_MAL=1.17e-3;   
deno=1+SUCc/KDCC_SUC+SUCm/KDCC_SUC+Pic/KDCC_PI+Pim/KDCC_PI...
    +MALc/KDCC_MAL+MALm/KDCC_MAL...
    +MALc*Pim/KDCC_MAL/KDCC_PI+MALm*Pic/KDCC_MAL/KDCC_PI...
    +SUCc*Pim/KDCC_SUC/KDCC_PI+SUCm*Pic/KDCC_SUC/KDCC_PI;
T_SUC_Pi=1*Tmax(5)/KDCC_SUC/KDCC_PI*(SUCc*Pim-Pic*SUCm)/deno;
%-------
%Pi-MAL   Pim+MALe---Pie+MALme  Antipoter 
deno=1+MALc/KDCC_MAL+MALm/KDCC_MAL+Pic/KDCC_PI+Pim/KDCC_PI...
    +SUCc/KDCC_SUC+SUCm/KDCC_SUC...
    +SUCc*Pim/KDCC_SUC/KDCC_PI+SUCm*Pic/KDCC_SUC/KDCC_PI...
    +MALc*Pim/KDCC_MAL/KDCC_PI+MALm*Pic/KDCC_MAL/KDCC_PI;

T_MAL_Pi=Tmax(6)/KDCC_MAL/KDCC_PI*(MALc*Pim-Pic*MALm)/deno;
%-------------------------------------------------------
%aKG-MAL   aKGm +MALe-aKGe+MALm  Antipoter 
K_aKG=0.24e-3;    
K_MAL=1e-3; 
deno=1+MALc/K_MAL+MALm/K_MAL+aKGm/K_aKG+aKGc/K_aKG...
    +aKGm*MALc/K_aKG/K_MAL+aKGc*MALm/K_aKG/K_MAL;
T_MAL_aKG=Tmax(7)/K_aKG/K_MAL*(aKGm*MALc-aKGc*MALm)/deno;
%-----------------------------------------------------------------
%%Carrier (TCC)   MALe+Hm+CITm-->He+CITe+MALm  Tricarboxylate 
%%
KTCC_MAL=0.25e-3;  
KTCC_CIT=1e-3;        
KTCC_H=1e-7;   
deno=1+MALc/K_MAL+MALm/K_MAL+CITm*Hm/KTCC_CIT/KTCC_H+CITc*Hc/KTCC_CIT/KTCC_H...
    +Hm*CITm*MALc/KTCC_CIT/KTCC_MAL/KTCC_H+Hc*CITc*MALm/KTCC_CIT/KTCC_MAL/KTCC_H;
T_MAL_HCIT=Tmax(8)/KTCC_MAL/KTCC_CIT/KTCC_H*(MALc*Hm*CITm-MALm*CITc*Hc)/deno;
%-------------------------------------
%PYR-H co-transporter between c and m (PYRH)
KPYRH_PYR=0.24e-3;      
KPYRH_H=1e-7;
deno=1+PYRm/KPYRH_PYR+PYRc/KPYRH_PYR+Hm/KPYRH_H+Hc/KPYRH_H+PYRm*Hm/KPYRH_H/KPYRH_PYR+PYRc*Hc/KPYRH_H/KPYRH_PYR;
T_PYRH=Tmax(9)/KPYRH_PYR/KPYRH_H*(PYRc*Hc-PYRm*Hm)/deno;
%----------------------------------------------------------
%Phosphate-H cotransporter (PIC)
K_PIC=9.4e-3;
K_PICH=1e-7;
deno=(1+Pic/K_PIC+Pim/K_PIC)*(1+Hm/K_PICH+Hc/K_PICH);
T_PIC=Tmax(10)/K_PICH/K_PIC*(Pic*Hc-Pim*Hm)/deno;
%--------------------------------------------------------------
% ANT
K_ATP=10e-6;
K_ADP=10e-6;
beta_ANT=0.6;
deno=(K_ATP*K_ADP+ADPc+ATPc*exp((beta_ANT-1)*F_con*dPsim/RT))*(K_ATP*K_ADP+ADPm+ATPm*exp(beta_ANT*F_con*dPsim/RT));
T_ANT =Tmax(11)*(exp(beta_ANT*dPsim*F_con/RT)*ADPc*ATPm-exp((beta_ANT-1)*dPsim*F_con/RT)*ADPm*ATPc)/deno;
%----------------------------------------------------
%GLUe*He-GLUm*Hm
KGLU_GLU=1e-3;
KGLU_H=1e-7;
deno=1+GLUc/KGLU_GLU+GLUm/KGLU_GLU+Hc/KGLU_H+Hm/KGLU_H+GLUc*Hc/KGLU_H/KGLU_GLU+GLUm*Hm/KGLU_H/KGLU_GLU;
T_GLUH=Tmax(12)/KGLU_H/KGLU_GLU*(GLUc*Hc-GLUm*Hm)/deno;

%------------------------------------------------
KAG_GLU=0.25e-3;       
KAG_ASP=0.12e-3;
KAG_H=1e-7;
%ASPc +HGLUm----ASPm +HGLUe
deno=1+Hm*GLUm/KAG_GLU/KAG_H+Hc*GLUc/KAG_GLU/KAG_H+ASPm/KAG_ASP+ASPc/KAG_ASP...
    +Hm*GLUm*ASPc/KAG_GLU/KAG_H/KAG_ASP+Hc*GLUc*ASPm/KAG_GLU/KAG_H/KAG_ASP;
T_ASP_GLU= Tmax(13)/KAG_GLU/KAG_ASP/KAG_H*(exp(dPsim.*F_con./(2*RT))*GLUm*Hm*ASPc-exp(-dPsim.*F_con./(2*RT))*ASPm*Hc*GLUc)/deno;
K_Hleak=1e-7;
if dPsim>1e-9
    T_HLEAK=Tmax(14)/K_Hleak*(Hc*exp(dPsim.*F_con./(2*RT))-Hm*exp(-dPsim.*F_con./(2*RT)));
else
    T_HLEAK=0;
end
%----------MA shuttle--------------------
Keq_APP_NADH=10^(pH_m-pH_c)*exp(dPsim.*F_con./RT);
deno=1+NADHc*NADm/1e-6+NADHm*NADc/1e-6;
   T_NADH=Tmax(15)*(Keq_APP_NADH*NADHc*NADm/1e-6-NADHm*NADc/1e-6)/deno;
%% Metabolic reaction and transport fluxes
%------------------------------------
rtfluxes0=[
J_HEX;
J_PGI;
J_PFK;
J_G3PF;
J_G3PD;
J_PHK;
J_PK;
J_LD;
J_PPP1;
J_PPP2;
J_GSR;
J_GP;
J_ATPASE;
J_AK;
J_MD2;
J_GOT2;
J_PDH;
J_CITS;
J_CITD;
J_AKGD;
J_SCAS;
J_NDK;
J_SUCD;
J_FUM;
J_MALD;
J_GOT;
J_CI;
J_CII;
J_CIII;
J_CIV;
J_CV;
T_GLCT_bc;
T_PYRT_bc;
T_LACT_bc;
T_PIT_bc;
T_SUC_Pi;
T_MAL_Pi;
T_MAL_aKG;
T_MAL_HCIT;
T_PYRH;
T_PIC;
T_ANT;
T_GLUH;
T_ASP_GLU;
T_HLEAK;
T_NADH;
J_ALA;];
rtfluxes=rtfluxes0;