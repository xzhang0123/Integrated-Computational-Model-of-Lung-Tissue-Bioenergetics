function obj=Lung_Parameter_Esimation_obj(x)
format long
global Flow F_con R_con Tem Vr Vb Vc Vm Vi
%Concentration unit is M
%%  Parameter Setup
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
    
Vr  =   55e-3;      %mL
Vb  =   0.66e-3; %mL
Vcell  =  0.67e-3; %mL  

Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
 Tem=310.15; %K      37 oC
 %Tem=298.15; %K      25 oC   25 degree
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
%------------------------------------------
%     iCytCoxi iCytCred iHi...
%     idPsim idPsip iO2 iR123e iR123m
Lung_dw=0.227; %each lung has 0.227g dry weight (1.33/5.87)
exp_glu05=[2.36 5.35 8.55 12.6 15.4]*Lung_dw;
exp_glu15=[4.67 11 19.5 27.7 36.4]*Lung_dw;
exp_glu3=[7.2 16.3 28.8 41.9 53.5]*Lung_dw;
exp_glu56=[10.8 23.5 35.7 56.2 69.5]*Lung_dw;
exp_glu10=[13.5 24.5 46.2 61.4 83.2]*Lung_dw;

%---------dGLC------lactate data------------
exp_lac05_exo=[13.1 24.8 30.4 32.6 33]*Lung_dw*0.18;  %Subtract endogenase part
exp_lac15_exo=[20.5 32.6 44.7 52.9 61]*Lung_dw*0.36;
exp_lac3_exo=[20.6 33.5 45.7 57.2 68.1]*Lung_dw*0.41;
exp_lac56_exo=[20 38.3 54.2 70.7 84.8]*Lung_dw*0.46;
exp_lac10_exo=[22.6 48.8 69 88.7 111]*Lung_dw*0.75;

exp_lac05=[13.1 24.8 30.4 32.6 33]*Lung_dw;  %Total
exp_lac15=[20.5 32.6 44.7 52.9 61]*Lung_dw; 
exp_lac3=[20.6 33.5 45.7 57.2 68.1]*Lung_dw;
exp_lac56=[20 38.3 54.2 70.7 84.8]*Lung_dw;
exp_lac10=[22.6 48.8 69 88.7 111]*Lung_dw;

exp_lac05_endo=exp_lac05-exp_lac05_exo;  %Subtract endogenase part
exp_lac15_endo=exp_lac15-exp_lac15_exo; %
exp_lac3_endo=exp_lac3-exp_lac3_exo;
exp_lac56_endo=exp_lac56-exp_lac56_exo;
exp_lac10_endo=exp_lac10-exp_lac10_exo;
%-----------dLAC data
factor_dLAC=1;
exp_dLAC0=[16 32.7 53.4 75 95.4]*Lung_dw/factor_dLAC;
%exp_dLAC05=[16 32.7 53.4 75 95.4]*Lung_dw/factor_dLAC*1.1;
exp_dLAC05=[16 32.7 53.4 75 95.4]*Lung_dw/factor_dLAC;
exp_dLAC2=[7.41 15.8 26.9 37 49.3]*Lung_dw/factor_dLAC;
exp_dLAC5=[4.77 12.5 19.9 27.7 35.8]*Lung_dw/factor_dLAC;
exp_dLAC20=[3.11 9.54 16.6 24.4 31.5]*Lung_dw/factor_dLAC;
%volumes
 %--------Determine paras to be estimated--------------------
 %Para(1:31) Vmaxs
  %Para(32:47) Tmaxs
  %Para(48) Mratio
  %Para(49) DH_factor
 Para=ones(1,51);
 
 %Vmax 3:8
  Para(3:8)=x(1:6);  %Vmax 3-8
Para(13)=x(7);  %ATPase
Para(46)=x(8);   %MA
Para(48)=x(9);     %Mass of Mito
Para(32:34)=x(10:12);  %Tmax
Para(49)=x(13);      % Vmax ALA
 Para(50)=x(14);      %Kctrl GLUT
% Para(51)=x(15);      %Kctrl ALA
Flow=12e-3;

% 
%  Para(1:9)=x(1:9);  %Vmax
% Para(32:34)=x(10:12);  %Tmax 1-3


options = odeset('NonNegative',[1:71]);
%------------------------------------------
  %%  Define t_step and t_final
t_step      =   20;   %min
%% Run Simulation
GLC_STEP=[0.5e-3 1.5e-3 3e-3 5.6e-3 10e-3 ];
LAC_STEP=[0e-3 0.5e-3 2e-3 5e-3 20e-3];
 time0=20;
 time1=100;
GLC_Mass=zeros(4,5);
LAC_Mass=zeros(4,5);
  IC=Set_Initial_Concentrations;
  %------------------------dGLC-Simulation------------------------
 for i=1:1:length(GLC_STEP)
      IC000=IC; 
    IC000(iGLCr)=GLC_STEP(i);
 [~,C00] = ode15s(@odeq,[0:0.1:time0],IC000,options,Para);
 IC1=C00(end,:);
[~,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
    if length(C(:,2))~=6
        obj=1e6;
        return
    end
%---------------------------
Tfluxes=zeros(length(C(:,1)),16)';
Rfluxes=zeros(length(C(:,1)),31)';
for istep=1:1:(length(C(1:end,1)))
    RTfluxes(:,istep)=fluxes(C(istep,:),Para);
end
Rfluxes=RTfluxes(1:31,:);
Tfluxes=RTfluxes(32:47,:);
%------------------------------------
GLC_Mass(i,:)=1e6*Vr*(IC1(iGLCr)-C((2:end),iGLCr)); %change unit from mol to umol
LAC_Mass(i,:)=1e6*Vr*(C((2:end),iLACr)-IC1(iLACr)); %change unit from mol to umol
% size(Rfluxes(7,:))
LAC_Mass_exo(i,:)=Rfluxes(7,2:end)./(Tfluxes(16,2:end)+Rfluxes(7,2:end)).*LAC_Mass(i,:);
LAC_Mass_endo(i,:)=Tfluxes(16,2:end)./(Tfluxes(16,2:end)+Rfluxes(7,2:end)).*LAC_Mass(i,:);
            if i==4

            for istep=1:1:(length(C(:,1)))
                RTfluxes(:,istep)=fluxes(C(istep,:),Para);
            end

  OCR=0.5*1e6*Rfluxes(30,3:end); %umol oxygen consumption rate, O2=0.5 CIV
% obj_ppp= 100*sum((PPP1-1/6*0.11*HEX).^2)/length(HEX); %11.9%

 obj_OCR= 200*sum((OCR-0.488).^2)./mean(OCR) ; %0.488 umol
    ATP_Total=C(3:end,iATPc)*Vc+C(3:end,iATPm)*Vm;
    ADP_Total=C(3:end,iADPc)*Vc+C(3:end,iADPm)*Vm;

 obj_ATPADP= 10*sum((ATP_Total./ADP_Total-4.831).^2)/length(ATP_Total)/(4.831^2); %  
%  obj_ATPAMP= sum((ATP_Total./AMP_Total-18.248).^2)/length(ATP_Total)/(18.248^2); %  
             end
 end
%  %GLC 5.6 mM is normal condition, so give these data more weights. *3
EXP_dGLC_GLC     =[exp_glu05 exp_glu15 exp_glu3 3*exp_glu56 exp_glu10];
EXP_dGLC_LAC     =[exp_lac05 exp_lac15 exp_lac3 3*exp_lac56 exp_lac10];
EXP_dGLC_LAC_ENDO=[exp_lac05_endo exp_lac15_endo exp_lac3_endo 3*exp_lac56_endo exp_lac10_endo];
EXP_dGLC_LAC_EXO =[exp_lac05_exo exp_lac15_exo exp_lac3_exo 3*exp_lac56_exo exp_lac10_exo];
%% %EXP_DATA_PYR=[0.1*exp_lac56];
 SIM_GLC=        [GLC_Mass(1,:) GLC_Mass(2,:) GLC_Mass(3,:) 3*GLC_Mass(4,:) GLC_Mass(5,:)];
 SIM_LAC=        [LAC_Mass(1,:) LAC_Mass(2,:) LAC_Mass(3,:) 3*LAC_Mass(4,:) LAC_Mass(5,:)];
 SIM_LAC_ENDO=   [LAC_Mass_endo(1,:) LAC_Mass_endo(2,:) LAC_Mass_endo(3,:) 3*LAC_Mass_endo(4,:) LAC_Mass_endo(5,:)];
SIM_LAC_EXO=     [LAC_Mass_exo(1,:) LAC_Mass_exo(2,:) LAC_Mass_exo(3,:) 3*LAC_Mass_exo(4,:) LAC_Mass_exo(5,:)];
obj_dGLC_GLC=sum(((SIM_GLC-EXP_dGLC_GLC)).^2)/length(SIM_GLC);
obj_dGLC_LAC=sum(((SIM_LAC-EXP_dGLC_LAC)).^2)/length(SIM_LAC);
obj_dGLC_LAC_EXO=sum(((SIM_LAC_EXO-EXP_dGLC_LAC_EXO)).^2)/length(SIM_LAC_EXO);
obj_dGLC_LAC_ENDO=sum(((SIM_LAC_ENDO-EXP_dGLC_LAC_ENDO)).^2)/length(SIM_LAC_ENDO);
%obj3=sum(((14*PYR_C-LAC_C)).^2)/length(PYR_C);
%% 
%---------dLAC simulations------------------------------------------;
Flow=14e-3;
 for i=1:1:length(LAC_STEP)
      IC000=IC; 
    IC000(iLACr)=LAC_STEP(i);
    IC000(iGLCr)=10e-3;
 [~,C00] = ode15s(@odeq,[0:0.1:20],IC000,options,Para);
 IC1=C00(end,:);
[~,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
    if length(C(:,2))~=6
        obj=1e6;
        return
    end
GLC_Mass(i,:)=1e6*Vr*(IC1(iGLCr)-C((2:end),iGLCr)); %change unit from mol to umol
LAC_Mass(i,:)=1e6*Vr*(C((2:end),iLACr)-IC1(iLACr)); %change unit from mol to umol
 end
 SIM_dLAC=[GLC_Mass(1,:) GLC_Mass(2,:) GLC_Mass(3,:) GLC_Mass(4,:) GLC_Mass(5,:)];
 exp_data_dLAC=[exp_dLAC0 exp_dLAC05 exp_dLAC2 exp_dLAC5 exp_dLAC20];
 obj_dLAC=10*sum(((SIM_dLAC-exp_data_dLAC)).^2)/length(SIM_dLAC);
 %% 
 %---------------------------------------------
 %dPO2 simulations;

  FEXP_PO2= [0.95 0.2 0.05 0.01 0.001 0.00006 ]*760;  %
FEXP_LAC2=[69.4 67.4 85.3 81.4 110 131 ]/60*0.227; %change unit to umol/min/lung
   FEXP_Ratio2=[8.7 9.3 10.8 12.3 16.8 20 ];
   FEXP_PYR2=[7.9 7.4 7.9 6.7 6.6 6.6]/60*0.227;
 dC_O2_step=0.003*FEXP_PO2/(22.4*1e6); %change Unit to: M

 time0=20;
 time1=40;
 for i=1:1:6
    
[~,C00] = ode15s(@odeq,[0:0.1:time0],IC,options,Para);
IC1=C00(end,:);
IC1(iO2)=dC_O2_step(i);
IC1(iGLCr)=10e-3;
[~,C1] = ode15s(@odeq,[0:10:time1],IC1,options,Para);

Tfluxes=zeros(length(C1(:,1)),16)';
for istep=1:1:(length(C1(1:end,1)))
    RTfluxes(:,istep)=fluxes(C1(istep,:),Para);
end
% Rfluxes=RTfluxes(1:31,:);
Tfluxes=RTfluxes(32:47,:);
T_LAC(i)=1e6*abs(Tfluxes(3,end)); %unit: umol/min  absolute volue
T_PYR(i)=1e6*abs(Tfluxes(2,end)); %unit: umol/min  absolute volue
Ratio_LP(i)=T_LAC(i)./T_PYR(i);
 end
  obj_dPO2_LAC=sum(((T_LAC-FEXP_LAC2)).^2)./mean(FEXP_LAC2);
 obj_dPO2_Ratio=sum(((Ratio_LP-FEXP_Ratio2)).^2)./mean(FEXP_Ratio2);
 obj_dPO2_PYR=20*sum(((T_PYR-FEXP_PYR2)).^2)./mean(FEXP_PYR2);

obj=obj_dLAC+obj_dGLC_GLC+10*obj_dGLC_LAC+obj_dPO2_LAC+obj_dPO2_Ratio...
    +obj_dPO2_PYR+obj_dGLC_LAC_EXO+obj_dGLC_LAC_ENDO+obj_OCR+obj_ATPADP
