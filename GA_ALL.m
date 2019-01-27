clear all
clc

tic

load General_Model_Parameters
toc
Algorithm=input('\n Enter 1 to use fmincon;\n Enter 2 to use Genetic algorithm;\n')


x0=dlmread('EstPara_final7.txt');
%x0=ones(1,15);
%------------------------
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end
%-----------------------------
LB=x0*0.1; % low bound
UB=x0*10;   %upper bound
% %--------------------------------
LB(9)=x0(9)*0.8;
UB(9)=x0(9)*1.2; % Mass of mito
if Algorithm==2

gaoptions = gaoptimset('MutationFcn',@mutationadaptfeasible);
gaoptions = gaoptimset(gaoptions,'PopulationSize',50,'Generations',100);
gaoptions = gaoptimset(gaoptions,'InitialPopulation',x0,'UseParallel',true);
gaoptions = gaoptimset(gaoptions,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr}, ...
    'Display','iter');
[Para_result,fval] = ga(@model_obj_GLU_PYR_LAC,length(LB),[],[],[],[],LB,UB,[],gaoptions)
elseif Algorithm==1
  options = optimoptions('fmincon','UseParallel',true,'Display','iter','Algorithm','interior-point');

Para_result = fmincon(@model_obj_GLU_PYR_LAC,x0,[],[],[],[],LB,UB,[],options)
end

dlmwrite('EstPara_final8.txt', Para_result,'precision','%10.14d');
