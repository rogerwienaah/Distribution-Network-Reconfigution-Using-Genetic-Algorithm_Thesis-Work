%***************************** Roger Wienaah, Randy Kofi Ansah, Cornelius Appiah *******************************
% ********************** Email: rogerwienaah@gmail.com, randyansah97@gmail.com, corneliusappiah39@gmail.com  **********************

clc
clear
clf 
close ALL

format short
global data1 NetPLoss st0 s t R X Z1 Vol_0
global Lweights Rweights tie tieSW NetVol %Vol_1
global Rweights1
global Lweights1



% --- tie switches -- Total num = 5
tie=[[8 20];[9 15];[12 22];[18 33];[25 29] ];  %-- branch nodes for tie sw
tieSW=[33; 34; 35;36;37];           % -- tie switch num
  [m, n]=size(tie);
  

  %---- s - source nodes for spanning tree
  %--- t - bus num
  s=[1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16 17 2  19 20 21 3  23 24 6 ...
      26 27 28 29 30 31 32];
  t=[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 ...
      27 28 29 30 31 32 33];
  s_t=[s' t'];
  
  % --- R = resistances and X = reactances of branch lines
  R =[0.0922 0.493  0.366  0.3811 0.819  0.1872 0.7114 1.03   1.044 ...
      0.1966 0.3744 1.468 0.5416  0.591  0.7463 1.289  0.732  0.164 ...
      1.5042 0.4095 0.7089 0.4512 0.898  0.896  0.203 0.2842  1.059 ...
      0.8042 0.5075 0.9744 0.3105 0.341];
     
  X =[0.047  0.2511 0.1864 0.1941 0.707   0.6188  0.2351  0.74    0.74 ...
      0.065  0.1238  1.155  0.7129 0.526  0.545  1.721  0.574  0.1565 ...
      1.3554  0.4784  0.9373  0.3083 0.7091 0.7011 0.1034  0.1447 ...
      0.9337 0.7006 0.2585 0.963  0.3619  0.5302];
  
  % ---- Lweights - real power of loads
  % --- original data from IEEE
  %--- rename to Lweights1 to use shuffle data
  Lweights   =[0 100 90 120 60 60 200 200 60 60 45 60 60 120 60 60 60 90 ...
            90 90 90 90 90 420 420 60 60 60 120 200 150 210 60];
 
  
  %--- shuffled load data      
  %LW=Lweights1(2:length(Lweights1));   
  %Lweights = [0 LW(randperm(length(LW)))];
  

  % --- Rweights - reactive power of load
  % --- original data from IEEE
  %--- rename to Rweights1 to use shuffle data
  Rweights  =[0 60 40  80 30 20  100 100 20 20 30 35 35 80 10 20 20 ...
           40 40 40 40 40 50 200 200 25 25 20 70 600 70 100 40];
  
  % -------This is use to insert DGs into the bus system
  %DG_pow()
  
  %--- shuffled load data
  %RW=Rweights1(2:length(Rweights1));   
  %Rweights = [0 RW(randperm(length(RW)))];
  
  
  %******************IEEE 33 Data for power flow Analysis*****************       
   data1=[s' t' R' X' Lweights(2:33)' Rweights(2:33)']; 
    st0 = [s', t'];
    
    % --- vol_0 contains the calculated pu Voltages and phase angles
    % --- Voltage_fun_span1 - does the load flow analysis
   [Vol_0,pow_L1,pow_X1,Z1]=Voltage_fun_span1(data1);
   NetPLoss=sum(pow_L1); % Net Real Power loss of IEEE33 network - standardized val based on loads
   NetVol=sum(Z1(:,5));
   Vol_1=Vol_0(:,1);
  
%********************************GA Optimization***************************
FitnessFcn=@Voltage_fun_span_optim;
ConstraintFunction=[];% @constraint_Voltage;
numberOfVariables = 5;          % used to define the search space and the number of decisions variables -- solution to return to select the best one -genes
rng(1,'twister')
IntCon = 1:5;
%% options setting
options = gaoptimset('TolFun',1e-12,'Generations',10  ,'Vectorized',...
    'off','PopulationSize',200, 'TolCon',1e-12,'Display','off',...
    'OutputFcn',@gaoutfun_span);
[x2,fval]=ga(@Voltage_fun_span_optim,numberOfVariables,[],[],[],[], ...
    [0 0 0 0 0],[1 1 1 1 1],ConstraintFunction,IntCon,options); 









