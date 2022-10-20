
% ----Plotting of results --------------------------------
function [state,options,optchanged] = gaoutfun_span(options,state,flag)
 global R X Vol_1 pow_L1 tie 
 global Z1 g data1 tieSW
 global s t Rweights Lweights 
 [m n]=size(tie);
  
optchanged = false;
switch flag
    case 'init'
        
    case {'iter','interrupt'}
         
    case 'done'
     
% %********plot Optimized IEEE33 graph with tie cycle switches*************

  s_t = [data1(:,1), data1(:,2)];
  st1=[s_t;g]; %Add tie switch data to IEEE 33 data
  s11=st1(:,1);
  t11=st1(:,2);
  R=data1(:,3);
  X=data1(:,4);

  ZBrL0=Z1(:,3); %sqrt(Z1(:,3).^2+Z1(:,4).^2);
 W=round([ ZBrL0; rand(numel(g(:,1)),1)],2); % tie swt random
 
%*********************plot graph with tie cycle switches******************
G = graph(s11,t11);
G.Edges.Weight=W;
figure (2)
subplot(121)
p = plot(G,'EdgeLabel',G.Edges.Weight);

%********Plot Optim.Minimum Tree with Krushal criteria*********************

 [T,pred] = minspantree(G,'Method','sparse');
 highlight(p,T);
 rootedTree = digraph(pred(pred~=0),find(pred~=0),[]);
TE=rootedTree.Edges{:,1}; % change Table to array
s3=TE(:,1);
t3=TE(:,2); 
W=round(ZBrL0,2);
G=graph(s3,t3);
G.Edges.Weight=W;
figure (2)
subplot(122)
p = plot(G,'EdgeLabel',G.Edges.Weight);
%p = plot(G);
title('IEEE33 Network Graph with OPTIM. Tie SW')
O_E=table2array(G.Edges(:,1));

%******************IEEE 33 Data for power flow Analysis********************
   s_t=[s' t'];
   data1=[s' t' R X Lweights(2:33)' Rweights(2:33)']; 
  
  [Vol_0,pow_L0,pow_X0,Z]=Voltage_fun_span1(data1);
  
%**********************Plot IEEE33 Power flow Results**********************
  figure (1)
  subplot(211),plot(1:numel(Vol_1(:,1)),Vol_1(:,1),'b',...
      1:numel(Vol_0(:,1)),Vol_0(:,1),'--')
  legend( 'Optimized', 'Base Case','location', 'southwest')
   title('IEEE33 Voltage Profile in p.u.')
   xlabel('Node')
   ylabel('p.u.')
  grid
  subplot(212),plot(1:numel(pow_L1(:,1)),pow_L1(:,1),'b', ...
      1:numel(pow_L0(:,1)),pow_L0(:,1),'--')
  sum1=sum(pow_L1);
  sum0=sum(pow_L0);
  legend( 'Optimized', 'Base Case','location', 'northeast')
  %title('IEEE33 Br. Power Losses')
  title({'IEEE33 Br. Power Losses ';['Base Case Loss: ',num2str(sum0),'KW  '...
      '  Optim. Loss: ',num2str(sum1),'KW']})
   xlabel('Node')
   ylabel('Power Loss KW')
  grid  
  
 %*************Plot IEEE33 Standard Span Tree******************************
  weights= Z(:,3);%sqrt(Z(:,3).^2 + Z(:,4).^2); % Power loss weights
   G = graph(s,t);
   G.Edges.Weight=round(weights,2);
  figure (3)
  subplot(121)
  %p = plot(G,'EdgeLabel',G.Edges.Weight);
  p = plot(G);
  title('IEEE33 Network Spanning Tree')
  S_E=[s' t'];
  E_str1 = string([s', t']);
  E_str=string([E_str1;tie(1:m,1:n)]);
 edgesUnique = cellstr(unique(E_str, 'rows'));
 index = find(strcmp(edgesUnique(:,1), edgesUnique(:,2)));
 if index~=0
 edgesUnique(index,:)=[];
 end
 s2=E_str(:,1);
 t2=E_str(:,2);

  weights_L=[ weights; rand(numel(tie(:,1)),1)]; %Power L.wts + tie dummy
 
 %**************Searching for sectional OFF switches********************** 
  comp=~ismember(S_E,O_E,'rows');
  stnd=S_E(comp,1)';
  %**************Searching for Tie ON/OFF switches************************* 
 % tieSWOFF=find(ismember(S_E,O_E,'rows')==0)
  tieSWOpen=tieSW(find(ismember(tie,g,'rows')));
  tieSWClosed=tieSW(find(~ismember(tie,g,'rows')));
  
 %******************Plot IEEE33 Graph with All Tie switches**************** 
load IEEE33EdgeTable
 s_t=IEEE33EdgeTable;
 %s_t=sortrows(s_t1,2)
Edge_L = string(s_t(:,3));
EdgeTable = table([s_t(:,1) s_t(:,2)],Edge_L, ...
    'VariableNames',{'EndNodes' 'Edge_L'});
Nodes = string(1:33)';
NodeTable = table(Nodes);
G =   graph(EdgeTable,NodeTable);
 figure (3)
 subplot(122)
plot(G,'NodeLabel',G.Nodes.Nodes,'EdgeLabel',G.Edges.Edge_L)
title('IEEE33 Network Graph with ALL SWs')

 %************************Printing Results*********************************
fprintf('\n')
fprintf('Tie Switches Closed:'),fprintf('%4g',tieSWOpen'), fprintf('\n')
fprintf('Switches Open:'),fprintf('%4g',stnd),...
fprintf('%4g',tieSWClosed), fprintf('\n')
 
end


