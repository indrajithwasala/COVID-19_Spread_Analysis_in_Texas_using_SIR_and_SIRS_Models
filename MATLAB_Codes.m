% MATLAB codes

% Data filtering and parameter estimation

clear all;
Texas = readtable('TexasT.csv');

Day=Texas{:,1};
t=Texas{:,2};
C_Tot=Texas{:,3};
D_Tot=Texas{:,4};

figure
subplot(1,2,1)
plot(Day,C_Tot,'r')
title('Cumulative covid cases in Texas')
subplot(1,2,2)
plot(Day,D_Tot,'b')
title('New daily covid cases in Texas')

% Calculate m
regfit=fitlm(t(1:244),log(C_Tot(1:244)/3),'Intercept',false);
regfit;
plot(regfit)

%m=0.071526

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;

% SIR and SIRS model fitting

%beta=1/28=0.035714  recovery rete
%alpha= m + beta = 0.107240  transmission rate 
%xsi=0.0; recovered changing rate to susceptible(This is zero in SIR modle)


a=0.142954;  % transmission rate
b=0.071428;  % recovery rete
c=0; % recovered changing rate to susceptible (This is zero in SIR modle)
N=29145505; %Total population in Texas
p=0; % The population had the vacation before the infection
y0=[29145502-0*N;3;p*N]; % s(0)=N-3=29,145,505-3 , i(0)=3 , r(0)=N*p
h=1; % time interval is one day

f=@(t,yy)[-(a/N)*yy(1)*yy(2)+c*yy(3);(a/N)*yy(1)*yy(2)-b*yy(2);
    b*yy(2)-c*yy(3)];
[T,yy]=ode45(f,0:h:600,y0); % Run 500 days

hold on
r1=plot(T,yy(:,1),'b');
r2=plot(T,yy(:,2),'r');
r3=plot(T,yy(:,3),'g');
hold off

xlabel({'Days','{(t)}'})
ylabel('Population ratios')
legend('Susceptible','Infected','Recovered')

%For ODE45 Infected Vs Susceptible plot
plot(yy(:,1),yy(:,2))
xlabel({'Susceptible ','{S(t)}'})
ylabel({'Infected','{I(t)}'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% Modeling school closure

a=0.107240;  % transmission rate
b=1/28; % recovery rete
c=0; % recovered changing rate to susceptible (This is zero in SIR modle)
N=29145505; %Total population in Texas
y0=[29145502-0*N;3;0*N]; % s(0)=N-3=29,145,505-3 , i(0)=3 , r(0)=N*p
h=1; % time interval is one day

f=@(t,yy)[-(a/N)*yy(1)*yy(2)+c*yy(3);(a/N)*yy(1)*yy(2)-b*yy(2);
    b*yy(2)-c*yy(3)];
f1=@(t,yy1)[0*yy1(1)*yy1(2)+c*yy1(3);0*yy1(1)*yy1(2)-b*yy1(2);
    b*yy1(2)-c*yy1(3)];
f2=@(t,yy2)[-(a/N)*yy2(1)*yy2(2)+c*yy2(3);(a/N)*yy2(1)*yy2(2)-b*yy2(2);
    b*yy2(2)-c*yy2(3)];

[T,yy]=ode45(f,0:h:249,y0); % Run 500 days
[T1,yy1]=ode45(f1,249:h:263,yy(250,:)); % Run 500 days
[T2,yy2]=ode45(f2,263:h:600,yy1(15,:));
yynew=[yy;yy1;yy2];
Tnew=[T;T1;T2];

hold on
r1=plot(Tnew,yynew(:,1),'b');
r2=plot(Tnew,yynew(:,2),'r');
r3=plot(Tnew,yynew(:,3),'g');
hold off

xlabel({'Days','{(t)}'})
ylabel('Population ratios')
legend('Susceptible','Infected','Recovered')
