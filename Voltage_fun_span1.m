% ----  Backward forward sweep -----
function [vbp,Plosskw,Qlosskw,Z0] = Voltage_fun_span1(xdata)

%MVAb = 1000;
%div = 1000;


%---loading data from IEEE ---
nf=xdata(:,1);
nt=xdata(:,2);
R_1=xdata(:,3);
Xl_1=xdata(:,4);
L_Act=[0; xdata(:,5)]; % --- actual power loss 
L_Reat=[0; xdata(:,6)]; % --- reactive power loss

br=length(R_1);
no=length(L_Act);

MVAb=100;%base MVA
KVb=11;% 11KV base
Zb=(KVb^2)/MVAb;% base impedance


% ---converting R_l and X_l to pu
for i=1:br
    R(i,1)= R_1(i)/Zb;
    Xl(i,1)=Xl_1(i)/Zb;
end

% --- converting P loss and Q loss to pu
for i=1:no
    P(i,1)=L_Act(i)/(1000*MVAb);
    Q(i,1)=L_Reat(i)/(1000*MVAb);
    %Qsh(i,1)=((loaddata(i,4))/(1000*MVAb));   
end

R;% line resistance
Xl;% line reactance
j=sqrt(-1);
Z=R+j*Xl;% line impedance
P;% load active power
Q;% load reactive power
Pt=sum(P);% total active power 
Qt=sum(Q);% total reactive power
%Pt
%Qt

C=zeros(br,no);% initilization for zero matrix (%---initializatin for the current in the lines)

for i=1:br
    a=nf(i);
    b=nt(i);
    for j=1:no
        if a==j
            C(i,j)=-1;
        end
        if b==j
            C(i,j)=1;
        end
    end
end
C;    
e=1;
for i=1:no
    d=0;
    for j=1:br
        if C(j,i)==-1
            d=1;
        end
    end
    if d==0
        endnode(e,1)=i;
        e=e+1;
    end
end
endnode;
h=length(endnode);
for j=1:h
    e=2;
    
    f=endnode(j,1);
   % while (f~=1)
   for s=1:no
     if (f~=1)
       k=1;  
       for i=1:br
           if ((C(i,f)==1)&&(k==1))
                f=i;
                k=2;
           end
       end
       k=1;
       for i=1:no
           if ((C(f,i)==-1)&&(k==1));
                f=i;
                g(j,e)=i;
                e=e+1;
                k=3;
           end            
       end
     end
   end
end
for i=1:h
    g(i,1)=endnode(i,1);
end
g;
w=length(g(1,:));
for i=1:h
    j=1;
    for k=1:no 
        for t=1:w
            if g(i,t)==k
                g(i,t)=g(i,j);
                g(i,j)=k;
                j=j+1;
             end
         end
    end
end
g;
for k=1:br
    e=1;
    for i=1:h
        for j=1:w-1
            if (g(i,j)==k) 
                if g(i,j+1)~=0
                    adjb(k,e)=g(i,j+1);            
                    e=e+1;
                else
                    adjb(k,1)=0;
                end
             end
        end
    end
end
adjb;
for i=1:br-1
    for j=h:-1:1
        for k=j:-1:2
            if adjb(i,j)==adjb(i,k-1)
                adjb(i,j)=0;
            end
        end
    end
end
adjb;
x=length(adjb(:,1));
ab=length(adjb(1,:));
for i=1:x
    for j=1:ab
        if adjb(i,j)==0 && j~=ab
            if adjb(i,j+1)~=0
                adjb(i,j)=adjb(i,j+1);
                adjb(i,j+1)=0;
            end
        end
        if adjb(i,j)~=0
            adjb(i,j)=adjb(i,j)-1;
        end
    end
end
adjb;
for i=1:x-1
    for j=1:ab
        adjcb(i,j)=adjb(i+1,j);
    end
end
b=length(adjcb);

% voltage current program  %--- backward / forward sweep part

for i=1:no
    vb(i,1)=1;
end
for s=1:10
for i=1:no
    nlc(i,1)=conj(complex(P(i,1),Q(i,1)))/(vb(i,1));
end
nlc;
for i=1:br
    Ibr(i,1)=nlc(i+1,1);
end
Ibr;
xy=length(adjcb(1,:));
for i=br-1:-1:1
    for k=1:xy
        if adjcb(i,k)~=0
            u=adjcb(i,k);
            %Ibr(i,1)=nlc(i+1,1)+Ibr(k,1);
            Ibr(i,1)=Ibr(i,1)+Ibr(u,1);
        end
    end      
end
Ibr;
for i=2:no
      g=0;
      for a=1:b 
          if xy>1
            if adjcb(a,2)==i-1 
                u=adjcb(a,1);
                vb(i,1)=((vb(u,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),Xl(i-1,1)))));
                g=1;
            end
            if adjcb(a,3)==i-1 
                u=adjcb(a,1);
                vb(i,1)=((vb(u,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),Xl(i-1,1)))));
                g=1;
            end
          end
        end
        if g==0
            vb(i,1)=((vb(i-1,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),Xl(i-1,1)))));
        end
end
s=s+1;
end
nlc;
Ibr;
vb;
vbp=[abs(vb) angle(vb)*180/pi];

for i=1:no
    va(i,2:3)=vbp(i,1:2);
end
for i=1:no
    va(i,1)=i;
end
va;


Ibrp=[abs(Ibr) angle(Ibr)*180/pi];
PL(1,1)=0;
QL(1,1)=0;

% Power losses ------------------------------------------objective function
for f=1:br
    Pl(f,1)=(Ibrp(f,1)^2)*R(f,1);
    Ql(f,1)=Xl(f,1)*(Ibrp(f,1)^2);
    PL(1,1)=PL(1,1)+Pl(f,1);
    QL(1,1)=QL(1,1)+Ql(f,1);
end

Plosskw=(Pl)*100000;
Qlosskw=(Ql)*100000;
PL=(PL)*100000;
QL=(QL)*100000;

voltage = vbp(:,1);
ang = vbp(:,2)*(pi/180);

for i=1:no
    if i==1
        Pg(i,1)=(Pt*100000)+PL;
        Qg(i,1)=(Qt*100000)+QL;
    else
        Pg(i,1)=0;
        Qg(i,1)=0;
    end
end

                       %display the power flow results
Vm=voltage';deltad=ang';Pgen=Pg';Qgen=Qg';Pd=P*100000';Qd=Q*100000';%Qinj=Qsh';

Pdt=sum(Pd);Qdt=sum(Qd);Pgent=sum(Pg);Qgent=sum(Qg);%Qinjt=sum(Qsh);


SLT = 0;
g=1;
for n = 1:no
busprt = 0;
   for L = 1:br;
       if busprt == 0

       busprt = 1;
       else
       end
       if nf(L)==n     k = nt(L);
       In = (Vm(n) - Vm(k))/Z(L);
       Ik = (Vm(k) - Vm(n))/Z(L);
       Snk = Vm(n)*conj(In)*100000;
       Skn = Vm(k)*conj(Ik)*100000;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       elseif nt(L)==n  k = nf(L);
       In = (Vm(n) - Vm(k))/Z(L);
       Ik = (Vm(k) - Vm(n))/Z(L);
       Snk = Vm(n)*conj(In)*100000;
       Skn = Vm(k)*conj(Ik)*100000;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       else
       end
         if nf(L)==n || nt(L)==n
         Knew=k;
          id(g,:)=[n Knew real(SL) imag(SL) Vm(n)];
         g=g+h;
         else
         end
  end
end
id(~any(id'),:)=[]; % Remove emply rows from Power data array
A=id; %[ActiveNLoad React.Nload ActiveLineLoss ReactiveLineLoss]   
[C,ia] = unique(A(:,4));
Z0=A(ia,:);

