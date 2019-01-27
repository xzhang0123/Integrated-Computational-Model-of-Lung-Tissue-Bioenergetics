function [EC,dflux]=Calculate_FCC(Para)
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

 
GLC_Mass=zeros(5,time1/t_step);
LAC_Mass=zeros(5,time1/t_step);

options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:71]);
 %---------dGLC---------------------------------------------
    IC=Set_Initial_Concentrations;
    IC(iGLCr)=5.6e-3;
[T00,C00] = ode15s(@odeq,[0:1:time0],IC,options,Para);
IC1=C00(end,:);
[~,C] = ode15s(@odeq,[0:t_step:time1],IC1,options,Para);
% GLC_Mass(i,:)=1e6*Vr*(C(1,iGLCr)-C(2:end,iGLCr)); %change unit from mol to umol
% LAC_Mass(i,:)=1e6*Vr*(C(2:end,iLACr)-C(1,iLACr)); %change unit from mol to umol
Mass_ATP=C(end,iATPm)*Vm+C(end,iATPc)*Vc;
Mass_ADP=C(end,iADPm)*Vm+C(end,iADPc)*Vc;
Mass_AMP=C(end,iAMPc)*Vc;
Tfluxes=zeros(length(C(:,1)),16)';

for istep=1:1:(length(C(1:end,1)))
    RTfluxes(:,istep)=fluxes(C(istep,:),Para);
end
Tfluxes=RTfluxes(32:47,:);
Rfluxes=RTfluxes(1:31,:);
 dflux=abs(Rfluxes(1,end))
 EC=(Mass_ATP+0.5*Mass_ADP)./(Mass_ATP+Mass_ADP+Mass_AMP);
%  ratio=Mass_ATP/Mass_ADP;
% size(Sim_ydata)
% size(Sim_Oxy_SUC)