function SSE=Calculate_SSE(Para,Edata)
format long
%Objective function for TCA cycle data at state 2
%%  Parameter Setup
global   Tem Flow F_con R_con Vr Vb Vc Vm Vi ROTi closed_system...
   iGLCr iPYRr iLACr iPir iHr...
   iGLCb iPYRb iLACb iPib iHb...
   iGLCc iG6Pc iF6Pc iF16BPc iGAPc iBPGc iPEPc iPYRc iLACc iPG6c...
   iR5Pc iMALc iOXAc iCITc iaKGc iSUCc iFUMc iGLUc iASPc iPic...
   iAMPc iADPc iATPc iNADHc iNADc iNADPHc iNADPc iGSSG iGSH iH2O2...
    iHc   iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iGLUm...
    iASPm iNADm iNADHm iACOAm iCOAm iUQm iUQH2m  iPim iADPm iATPm...
    iFADm iFADH2m iHm iCytCoxi iCytCred iHi idPsim idPsip iO2 iR123e iR123m

%% GLC and LAC simulation
t_step      =   20;   %min
time0=20; 
time1=100;
GLC_STEP=[0.5e-3 1.5e-3 3e-3 5.6e-3 10e-3 ];
LAC_STEP=[0e-3 0.5e-3 2e-3 5e-3 20e-3];
FEXP_PO2= [0.95 0.2 0.05 0.01 0.001 0.00006]*760;  %
 dC_O2_step=0.003*FEXP_PO2/(22.4*1e6); %change Unit to: M
 
GLC_Mass=zeros(5,time1/t_step);
LAC_Mass=zeros(5,time1/t_step);

options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);
 %---------dGLC---------------------------------------------
for i=1:1:length(GLC_STEP)
    IC=Set_Initial_Concentrations;
    IC(iGLCr)=GLC_STEP(i);
[T00,C00] = ode15s(@odeq,[0:1:time0],IC,options,Para);
IC1=C00(end,:);
[~,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
GLC_Mass(i,:)=1e6*Vr*(C(1,iGLCr)-C(2:end,iGLCr)); %change unit from mol to umol
LAC_Mass(i,:)=1e6*Vr*(C(2:end,iLACr)-C(1,iLACr)); %change unit from mol to umol
end
SIM_dGLC=[GLC_Mass(1,:) GLC_Mass(2,:) GLC_Mass(3,:) GLC_Mass(4,:) GLC_Mass(5,:) ...
    LAC_Mass(1,:) LAC_Mass(2,:) LAC_Mass(3,:) LAC_Mass(4,:) LAC_Mass(5,:) ];
%-----------------dLAC simulations------------------

 for i=1:1:length(LAC_STEP)
      IC000=IC; 
    IC000(iLACr)=LAC_STEP(i);
 [~,C00] = ode15s(@odeq,[0:0.1:20],IC000,options,Para);
 IC1=C00(end,:);
[~,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);

GLC_Mass(i,:)=1e6*Vr*(IC1(iGLCr)-C((2:end),iGLCr)); %change unit from mol to umol
LAC_Mass(i,:)=1e6*Vr*(C((2:end),iLACr)-IC1(iLACr)); %change unit from mol to umol
 end
 SIM_dLAC=[GLC_Mass(1,:) GLC_Mass(2,:) GLC_Mass(3,:) GLC_Mass(4,:) GLC_Mass(5,:)];
%--------------dPO2----------------------------------
 for i=1:1:6
    
[~,C00] = ode15s(@odeq,[0:0.1:time0],IC,options,Para);
IC1=C00(end,:);
IC1(iO2)=dC_O2_step(i);
[~,C1] = ode15s(@odeq,[0:t_step:60],IC1,options,Para);

Tfluxes=zeros(length(C1(:,1)),16)';

for istep=1:1:(length(C1(1:end,1)))
    RTfluxes(:,istep)=fluxes(C1(istep,:),Para);
end
Tfluxes=RTfluxes(32:47,:);
T_LAC(i)=1e6*abs(Tfluxes(3,end)); %unit: umol/min  absolute volue
T_PYR(i)=1e6*abs(Tfluxes(2,end)); %unit: umol/min  absolute volue
Ratio_LP(i)=T_LAC(i)./T_PYR(i);
 end
SIM_dPO2=[T_LAC T_PYR log10(Ratio_LP)]; %increase T_PYR weighting?
 SIM=[SIM_dGLC SIM_dLAC SIM_dPO2];
 SSE=sum((SIM-Edata).^2);
% size(Sim_ydata)
% size(Sim_Oxy_SUC)