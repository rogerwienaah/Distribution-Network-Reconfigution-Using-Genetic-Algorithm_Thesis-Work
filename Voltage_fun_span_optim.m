% Minimum Spanning Tree Algorithm
function err = Voltage_fun_span_optim(x)
 
global data1 NetPLoss s1 t1 Z1 u g tie Z2  ZBrL0 
global  s2 t2 R1 X1 Vol_1 pow_L1 pow_X1 data2 x0 
  
  x0=x';
  x2f=x0.*tie(:,:);
  g1=x2f(:,1);
  g2=x2f(:,2);
  g=[g1(g1~=0) g2(g2~=0)];
  u = g; % Otimized tie SW
  s_t = [data1(:,1), data1(:,2)];
  st1=[s_t;g]; %Add tie switch data to IEEE 33 data
  s1=st1(:,1);
  t1=st1(:,2);
  R=data1(:,3);
  X=data1(:,4);

  ZBrL0=Z1(:,3); %sqrt(Z1(:,3).^2+Z1(:,4).^2);
 W1=round([ ZBrL0; rand(numel(g(:,1)),1)],2); % tie swt random
%*********************plot graph with tie cycle switches******************
G = graph(s1,t1);
G.Edges.Weight=W1;
%*************Plot Minimum Tree with Krushal criteria********************** 
[T,pred] = minspantree(G,'Method','sparse');

rootedTree = digraph(pred(pred~=0),find(pred~=0),[]);

%*******************Minimum Tree Egdes data********************************
TE=rootedTree.Edges{:,1}; % change Table to array
s2=TE(:,1);
t2=TE(:,2);
netdata0=[s2, t2];
netdata0=sortrows(netdata0,2);

%********Inserting tie line data to Min. Span Tree Edge data***************
v=length(g(:,1));
for y=1:v
    f0=ismember(netdata0,g(y,:));
    f1=find(f0(:,1)&f0(:,2)==1);  
        R(f1,1)=2;
        X(f1,2)=2;
end

R1=R;
X1=X;

 %************************Network Power Loss and Voltage****************
 netWt=[data1(:,5) data1(:,6)];
 data2=[netdata0(:,1) netdata0(:,2) R1 X1 netWt(:,1) netWt(:,2)]; 
[Vol_1,pow_L1,pow_X1,Z2]=Voltage_fun_span1(data2);
err= (sum(Z2(:,3)-NetPLoss));%+sum(Vol_1(:,1)-15));%15
