
function  IC=Set_Initial_Concentrations
%%  Set default initial conditions (concentration, [M] )

%Dry weight of each isolated lung=1.33/5.87=0.227 g 
format long

  
Vcell  =  0.67e-3; %mL  
Vm=1/51.0714*Vcell;% mL  2% 
Vc=50/51.0714*Vcell;
Vi=0.0724/51.0714*Vcell;
%----------------calculate pool concentrations for NADH C,ANP, COA----------------
Vmatrix=1e-6;  %1uL/mg mitochondria
%Mito_protein_mass=Mito_density*Vm; % =3.57 each lung has 3.57mg mitochondria protein
Poolm_NADH=1.73e-9/Vmatrix;  %Concentration, Total NADH=1.73nmol/mg of mito protein. 
Poolm_CytC=0.33e-9/Vmatrix;
Poolm_ANP=6.4e-9/Vmatrix;   %/ Vi
Poolm_COA=0.925e-9/Vmatrix;
%-----------------------------------------------
GLCr=5.6e-3;
PYRr=0;
LACr=0;
Pir=2.5e-3;
Hr=10^(-7.4);
%------------
GLCb=GLCr;
PYRb=PYRr;
LACb=LACr;
Pib=Pir;
Hb=Hr;
%------------------
GLCc=0.56e-3;  %[1]
G6Pc=0.278e-3;  %[1]
F6Pc=0.0486e-3;
F16BPc=0.0667e-3;
GAPc=0.0889e-3;
BPGc=88.9e-6;
PEPc=5.39e-6;
PYRc=0.1e-3; % model paper
LACc=1e-3;   %ref= model paper
PG6c=0.3e-3;
R5Pc=0.01e-3;
MALc=2e-3;
OXAc=2e-3;
CITc=2e-3;
aKGc=2e-3;
SUCc=2e-3;
FUMc=2e-3;
GLUc=6e-3;
ASPc=6e-3;
Pic=2.9e-3;
%AMPc=90.8e-9/Vc;  %0.1384  %whole lung AMP content:90.8-138nmol/lung; another ref=energy status of the rat lung after exposure to ++Po2
AMPc=138e-9/Vc;
ADPc=(204.3-11.4)*1e-9/Vc;   %Ref=Akai et al. 204.3  =0.294 mM
ATPc=(1589-11.4)*1e-9/Vc;    %                   =   2.4mM
NADHc=0.68e-6;
NADc=340e-6;
NADPHc=22.7e-9/Vc;  %From Akai et al.
NADPc=15.9e-9/Vc;   
GSSG=136.2e-9/Vc;
GSH=2724e-9/Vc;
H2O2=0.03e-3;
Hc=10^(-7.1);
%-------------------------
PYRm =0.025e-3;           
OXAm = 0.001e-3;     
CITm =3e-3;        
aKGm = 3e-3;      
SCAm = 0.1*Poolm_COA;       
SUCm = 2e-3;        
FUMm=3e-3;
MALm = 3e-3;        
GLUm= 6e-3;
ASPm =6e-3;
NADm=0.25*Poolm_NADH;
NADHm=0.75*Poolm_NADH;
ACOAm = 0.4*Poolm_COA;    
COAm = 0.5*Poolm_COA;   
UQpool=1.35e-3;
UQm=0.5*UQpool;  
UQH2m=0.5*UQpool;     
Pim = 2e-3;        
ADPm = 0.5*Poolm_ANP;   
ATPm =0.5*Poolm_ANP;     
Hm = 10^(-7.6);   
FADm =0.6* 0.7e-3;       
FADH2m = 0.4*0.7e-3;    
 CytCoxi = 0.9*Poolm_CytC;  
 CytCred = 0.1*Poolm_CytC; 
Hi=10^(-7.1);
dPsim = 130;      %Unit: millivolts
dPsip=40;
O2 = 0.296e-3;      
R123e=200e-9;      
R123m=0;
%% save initialconcentrations
 initialconditions = [
GLCr;
PYRr;
LACr;
Pir;
Hr;
%--------------------
GLCb;
PYRb;
LACb;
Pib;
Hb;
%----------------------
GLCc;
G6Pc;
F6Pc;
F16BPc;
GAPc;
BPGc;
PEPc;
PYRc;
LACc;
PG6c;
R5Pc;
MALc;
OXAc;
CITc;
aKGc;
SUCc;
FUMc;
GLUc;
ASPc;
Pic;
AMPc;
ADPc;
ATPc;
NADHc;
NADc;
NADPHc;
NADPc;
GSSG;
GSH;
H2O2;
Hc;
%------------------
PYRm;
OXAm;
CITm;
aKGm;
SCAm;
SUCm;
FUMm;
MALm;
GLUm;
ASPm;
NADm;
NADHm;
ACOAm;
COAm;
UQm;
UQH2m;
Pim;
ADPm;
ATPm;
FADm;
FADH2m;    
Hm;
%------------------
CytCoxi;
CytCred;
Hi;
dPsim;
dPsip; %set to constant
O2;    %constant
R123e;  %not used for intact lung model
R123m;   %not usded for intact lung model
   ];
IC=initialconditions;
