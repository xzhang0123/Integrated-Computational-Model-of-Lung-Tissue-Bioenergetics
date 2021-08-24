function dxdt = odeq(t,IC,Para1)

global Flow R_con Tem Vr Vb Vc Vm Vi 
%   ode equations
format long
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
%% Load data from global variable       
CO2=1.32e-5; %CO2-Henry's Law: 3.87e-4 atm/(29.41 atm/M)
RT=R_con*Tem;
buffer_m= 2.303*10^(-7.6)/40000;
buffer_c=2.303*10^(-7.1)/2500;
buffer_i=2.303*10^(-7.1)/25000;
rtfluxes=fluxes(IC,Para1);
J_HEX=rtfluxes(1);
J_PGI=rtfluxes(2);
J_PFK=rtfluxes(3);
J_G3PF=rtfluxes(4);
J_G3PD=rtfluxes(5);
J_PHK=rtfluxes(6);
J_PK=rtfluxes(7);
J_LD=rtfluxes(8);
J_PPP1=rtfluxes(9);
J_PPP2=rtfluxes(10);
J_GSH1=rtfluxes(11);
J_GSH2=rtfluxes(12);
J_ATPASE=rtfluxes(13);
J_AK=rtfluxes(14);
J_MD2=rtfluxes(15);
J_GOT2=rtfluxes(16);
J_PDH=rtfluxes(17);
J_CITS=rtfluxes(18);
J_CITD=rtfluxes(19);
J_AKGD=rtfluxes(20);
J_SCAS=rtfluxes(21);
J_NDK=rtfluxes(22);
J_SUCDH=rtfluxes(23);
J_FUM=rtfluxes(24);
J_MALD=rtfluxes(25);
J_GOT=rtfluxes(26);
J_CI=rtfluxes(27);
J_CII=rtfluxes(28);% TAKE OUT CII
J_CIII=rtfluxes(29);
J_CIV=rtfluxes(30);
J_CV=rtfluxes(31);
T_GLCT_bc=rtfluxes(32);
T_PYRT_bc=rtfluxes(33);
T_LACT_bc=rtfluxes(34);
T_PIT_bc=rtfluxes(35);
T_SUC_Pi=rtfluxes(36);
T_MAL_Pi=rtfluxes(37);
T_MAL_aKG=rtfluxes(38);
T_MAL_HCIT=rtfluxes(39);
T_PYRH_cm=rtfluxes(40);
T_PIC=rtfluxes(41);
T_ANT=rtfluxes(42);
T_GLUH=rtfluxes(43);
T_ASP_GLU=rtfluxes(44);
T_HLEAK=rtfluxes(45);
T_NADH=rtfluxes(46);
T_ALA=rtfluxes(47);
%----------------------
GLCr=IC(iGLCr);
PYRr=IC(iPYRr);
LACr=IC(iLACr);
Pir=IC(iPir);
Hr=IC(iHr);
%--------------------
GLCb=IC(iGLCb);
PYRb=IC(iPYRb);
LACb=IC(iLACb);
Pib=IC(iPib);
Hb=IC(iHb);
%% Net Equations
dxdt(iGLCr)=(Flow/Vr)*(GLCb - GLCr);
dxdt(iPYRr)=(Flow/Vr)*(PYRb - PYRr);
dxdt(iLACr)= (Flow/Vr)*(LACb - LACr);
dxdt(iPir)= (Flow/Vr)*(Pib  - Pir);
dxdt(iHr)=0;
%recirculating: Concentrations are constant
% dxdt(iGLCr)=0;
% dxdt(iPYRr)=0;
% dxdt(iLACr)= 0;
% dxdt(iPir)=0;
%--------------------------------------
dxdt(iGLCb)= (Flow*(GLCr - GLCb) - T_GLCT_bc)/(Vb);
dxdt(iPYRb)=(Flow*(PYRr - PYRb) - T_PYRT_bc)/(Vb);
dxdt(iLACb)=(Flow*(LACr - LACb) - T_LACT_bc)/(Vb);
dxdt(iPib)=(Flow*(Pir - Pib) - T_PIT_bc )/(Vb);
dxdt(iHb)=0;
%--------------------------------------------
dxdt(iGLCc)=(T_GLCT_bc-J_HEX)/Vc;
dxdt(iG6Pc)=(J_HEX - J_PGI - (1/6)*J_PPP1)/Vc;
dxdt(iF6Pc)=(J_PGI - J_PFK)/Vc;
dxdt(iF16BPc)= ( J_PFK - J_G3PF)/Vc;
dxdt(iGAPc)= (2*J_G3PF-J_G3PD)/Vc;
dxdt(iBPGc)= (J_G3PD-J_PHK)/Vc;
dxdt(iPEPc)= (J_PHK-J_PK)/Vc;
dxdt(iPYRc)= (J_PK-J_LD+T_PYRT_bc-T_PYRH_cm+T_ALA)/Vc;
dxdt(iLACc)=(J_LD+T_LACT_bc)/Vc;
dxdt(iPG6c) = 0;
dxdt(iR5Pc) =0;
dxdt(iMALc) = 0*(-T_MAL_HCIT)/Vc;
dxdt(iOXAc) =(-J_MD2-J_GOT2)/Vc;
dxdt(iCITc) = (+T_MAL_HCIT)/Vc;
dxdt(iaKGc) = (+T_MAL_aKG+J_GOT2)/Vc;
dxdt(iSUCc) = 0*(-T_SUC_Pi)/Vc;
dxdt(iFUMc) = (0)/Vc;    %NO REACTION
dxdt(iGLUc) = (-T_GLUH+T_ASP_GLU-J_GOT2)/Vc;
dxdt(iASPc) = (-T_ASP_GLU+J_GOT2)/Vc;
dxdt(iPic)  =(T_PIT_bc+T_SUC_Pi +T_MAL_Pi-J_G3PD+J_ATPASE-T_PIC)/Vc;
dxdt(iAMPc) =-J_AK/Vc;
dxdt(iADPc) =(J_HEX+J_PFK-J_PHK-J_PK+2*J_AK-T_ANT+J_ATPASE)/Vc;
dxdt(iATPc)  =(-J_HEX-J_PFK+J_PHK+J_PK-J_AK+T_ANT-J_ATPASE)/Vc;
 dxdt(iNADHc) = (-T_NADH+J_G3PD-J_LD-J_MD2)/Vc;
dxdt(iNADc)   =-dxdt(iNADHc);
dxdt(iNADPHc)=(+2*J_PPP1-J_GSH1)/Vc;
dxdt(iNADPc)=-dxdt(iNADPHc);
dxdt(iGSSG)=(-J_GSH1+J_GSH2)/Vc;
dxdt(iGSH) =(2*J_GSH1-2*J_GSH2)/Vc;
dxdt(iH2O2)=0;
dxdt(iHc) =0;
%Mitochondria region
dxdt(iPYRm) = (T_PYRH_cm - J_PDH)/Vm;
dxdt(iOXAm)= (-J_CITS + J_MALD+J_GOT)/Vm;
dxdt(iCITm) = (-T_MAL_HCIT + J_CITS - J_CITD)/Vm;
dxdt(iaKGm)  = (J_CITD - J_GOT - J_AKGD - T_MAL_aKG)/Vm;
dxdt(iSCAm)  = (J_AKGD - J_SCAS)/Vm;
dxdt(iSUCm)  = (+T_SUC_Pi + J_SCAS - J_SUCDH)/Vm;
dxdt(iFUMm)  =(+J_SUCDH-J_FUM)/Vm;
dxdt(iMALm)  = (J_FUM - J_MALD  +T_MAL_Pi + T_MAL_aKG + T_MAL_HCIT)/Vm;
dxdt(iGLUm) = (J_GOT - T_ASP_GLU + T_GLUH )/Vm;
dxdt(iASPm) = (-J_GOT+T_ASP_GLU)/Vm;
dxdt(iNADm) = (-T_NADH-J_PDH - J_CITD - J_AKGD - J_MALD + J_CI)/Vm;
dxdt(iNADHm) = -dxdt(iNADm);
dxdt(iACOAm)=(J_PDH-J_CITS)/Vm;
dxdt(iCOAm)=(- J_PDH+J_CITS + J_SCAS - J_AKGD )/Vm;
dxdt(iUQm)=(-J_CI-J_SUCDH+J_CIII)/Vm;
dxdt(iUQH2m)=-dxdt(iUQm);
dxdt(iPim)  = (-T_SUC_Pi - T_MAL_Pi - J_SCAS - J_CV +T_PIC  )/Vm;
dxdt(iADPm) = (T_ANT - J_SCAS - J_CV )/Vm;
dxdt(iATPm) = (J_CV + J_SCAS - T_ANT )/Vm;
dxdt(iFADm)  = 0*(J_CII-J_SUCDH)/Vm;   %TAKE OUT FAD
dxdt(iFADH2m)=0* -dxdt(iFADm);         %TAKE OUT
dxdt(iHm)   = (buffer_m*( -T_MAL_HCIT+T_PYRH_cm+ T_PIC+T_HLEAK+  - (4+1)*J_CI - 2*J_CIII-4*J_CIV+(3-1)*J_CV))/Vm;
dxdt(iCytCoxi)=(-2*J_CIII+2*J_CIV)/Vi;
dxdt(iCytCred)=-dxdt(iCytCoxi);
dxdt(iHi)   = buffer_i*(4*J_CI+4*J_CIII+2*J_CIV-3*J_CV-T_HLEAK-T_PIC)/Vi;
Cimm=6.75*1e-6*(Vm+Vi);       
dxdt(idPsim)=1/Cimm*(4*J_CI+2*J_CIII+4*J_CIV-3*J_CV-T_HLEAK-T_ANT-T_NADH);
dxdt(idPsip)=0;
dxdt(iO2)=0;%open system, oxygen is constant
dxdt(iR123e)=0;
dxdt(iR123m)=0;
dxdt=dxdt';

end

