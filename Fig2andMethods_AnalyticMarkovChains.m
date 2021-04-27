% Figure 6
% Analytic treatment of an ordinary differential equation representation
% of the futile cycle promoter
% 


% a colormap 
cmap = GetColorMap('hsv',8);

%% Fig 6B
% parameters 
k1 = .003; % Rate of conversion from state 0 to state 1. 
k2 = .5; % Rate of conversion from state 1 to state 2. 
tau1 = 20; % Rate of conversion from state 1 to state 0.
tau2 = .05; % Rate of conversion from state 2 to state 1.
hdac = 2; % Rate of tag removal. 
r0 = .01; % Basal tagging rate (state 0)
r1 = .5; % Stimulated tagging rate (state 1) 
r2 =5.0; % Maximally stimulated tagging rate (state 2)
Et= 20; % Maximum tag-transferases associated with promoter
%
syms a
loopRates = 0:.2:20;
nSteps = length(loopRates);
lows = zeros(nSteps,1);
highs = zeros(nSteps,1);
for i=1:nSteps
    dadt = (Et*(k1*k2*r2*a.^2 + k1*r1*tau2*a + r0*tau1*tau2))./(k1*k2*a.^2 + k1*tau2*a + tau1*tau2) +loopRates(i) -a*hdac;
    soln = vpasolve(dadt==0,a,[0,60]);
    if ~isempty(soln)
        lows(i) = min(soln);
        highs(i) = max(soln);
    end
end
%
f2 = figure(2); clf;
plot(loopRates,lows,'bo','MarkerSize',4); hold on;
plot(loopRates,highs,'r.');
plot(loopRates,lows,'b-'); hold on;
plot(loopRates,highs,'r-');
xlabel('contact frequency');
ylabel('transcription rate');
set(gcf,'color','w');
xlim([-1,19]); ylim([0,70]);

%% Fig 6C: tagConc_vs_tagRate for diff enh loop rates
f2 = figure(2); clf; 
c = 0;
for loopRate =linspace(0,14,8)
    c=c+1;
    p = p1 +loopRate;
    plot(a,p-hdac*a,'color',cmap(c,:)); hold on;
end
xlabel('tag concentration [a]');
ylabel('tagging rate [da/dt]');
title('loop Rate Scan');
legend(  cellstr(num2str(linspace(0,14,8)',3)));
plot([0,50],[0,0],'k-'); hold on;
set(gcf,'color','w');

%% Fig 6E: changing N

%---------   n=2 ---------------%
% states
% ac0 <=> ac1 <=> ac2  <=> ac3
% Pr <=> Pr-E1 <=> Pr-E2
SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Images\');
syms E0 E1 E2 Et tau1 tau2 k1 k2 r0 r1 r2 hdac a;
consE = Et -(E0 + E1 + E2);
dE0 = -k1*a*E0 + tau1*E1;
dE1 =  k1*a*E0 + tau2*E2 -k2*a*E1 -tau1*E1;
dE2 =  k2*a*E1 - tau2*E2;
nullE = solve([dE0==0,dE1==0,dE2==0,consE==0],[E0,E1,E2]);
nullE.E0
nullE.E1
nullE.E2
da = r0*E0 + r1*E1 + r2*E2 - hdac*a;

simplify( r0*nullE.E0 + r1*nullE.E1 + r2*nullE.E2)
n2 = (Et*(k1*k2*r2*a^2 + k1*r1*tau2*a + r0*tau1*tau2))/(k1*k2*a^2 + k1*tau2*a + tau1*tau2)
% n=1
SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Images\');
syms E0 E1 E2 Et tau1 tau2 k1 k2 r0 r1 r2 hdac a;
consE = Et -(E0 + E1);
dE0 = -k1*a*E0 + tau1*E1;
dE1 =  k1*a*E0  -tau1*E1;
nullE = solve([dE0==0,dE1==0,consE==0],[E0,E1]);
nullE.E0
nullE.E1
simplify( r0*nullE.E0 + r1*nullE.E1 )

n1 = (Et*(r0*tau1 + a*k1*r1))/(tau1 + a*k1)

% n=3
syms E0 E1 E2 E3 Et tau1 tau2 tau3 k1 k2 k3 r0 r1 r2 r3 hdac a;
consE = Et -(E0 + E1 + E2 + E3);
dE0 = -k1*a*E0 + tau1*E1;
dE1 =  k1*a*E0 + tau2*E2 -k2*a*E1 -tau1*E1;
dE2 =  k2*a*E1 + tau3*E3 -k3*a*E2 - tau2*E2;
dE3 =  k3*a*E2 - tau3*E3;
nullE = solve([dE0==0,dE1==0,dE2==0,dE3==0,consE==0],[E0,E1,E2,E3]);
da = r0*E0 + r1*E1 + r2*E2 + r3*E3 - hdac*a;

simplify( r0*nullE.E0 + r1*nullE.E1 + r2*nullE.E2+ r3*nullE.E3)
n3 = (Et*(k1*k2*k3*r3*a^3 + k1*k2*r2*tau3*a^2 + k1*r1*tau2*tau3*a + r0*tau1*tau2*tau3))/(k1*k2*k3*a^3 + k1*k2*tau3*a^2 + k1*tau2*tau3*a + tau1*tau2*tau3)

% Changing N, plot results
k3 = .005; tau3 =.5; 
k1 = .0022; % higher steepens the curve, pulling the Km to lower values. 
k2 = .5; % higher steepens the curve, pulling the Km to lower values
tau1 = 15; % just scales the ks
tau2 = .05;
hdac = 2; % just scales the rs. 
r0 = .01; % .7; % 3/1000/.005; % sets behavior at ac=0, effects loopRate contribution
r1 = .5; % not super sensitive 
r2 =5.0; % large values for low ks pushes Km to higher values
r3 = 4; 
loopRate = 4; % .65 ;% .95; % the linear shift to the sigmoid. % 20% change, 3.5 to 4.5  
Et= 20;
a = 0:0.1:50; 

n1 = (Et*(r0*tau1 + a*k1*r1))./(tau1 + a*k1);
n2 = (Et*(k1*k2*r2*a.^2 + k1*r1*tau2*a + r0*tau1*tau2))./(k1*k2*a.^2 + k1*tau2*a + tau1*tau2);
n3 = (Et*(k1*k2*k3*r3*a.^3 + k1*k2*r2*tau3*a.^2 + k1*r1*tau2*tau3*a + r0*tau1*tau2*tau3))./(k1*k2*k3*a.^3 + k1*k2*tau3*a.^2 + k1*tau2*tau3*a + tau1*tau2*tau3);

f2 = figure(2); clf;
plot(a,n1+loopRate -hdac*a); hold on;
plot(a,n2+loopRate -hdac*a);
plot(a,n3+loopRate -hdac*a);
ylim([-10,10]);
xlabel('tag concentration [a]');
ylabel('tagging rate [da/dt]');
title('n Scan');
legend(  '1','2','3');
plot([0,50],[0,0],'k-'); hold on;
set(gcf,'color','w');

%% Fig 6f: r0 scan
figure(3); clf;%  
f5 = figure(5); clf;
plot([0,50],[0,0],'k-'); hold on;
c = 0;
Et = 20;
r2 = 5.0;  r1 = .5; r0=.01; 
loopRate = 5; k1=.003;  k2 = .5;
tau1 = 20; tau2 = .05;
vs =linspace(.005,10,8);%  5*logspace(-2,0,8) % 
hds = fliplr(linspace(.02,2,8));
a = linspace(0,50);
for r0 =vs % 
    c=c+1;
    p = (Et*(k1*k2*r2./(tau1*tau2)*a.^2 + k1*r1*tau2./(tau1*tau2)*a + r0))...
        ./(k1*k2./(tau1*tau2)*a.^2 + k1*tau2./(tau1*tau2)*a + 1  )...
        + loopRate 
    p1 = p -hds(c)*a;
   figure(3); plot(a,p1,'color',cmap(c,:)); hold on;
   figure(5); plot(a,p,'color',cmap(c,:)); hold on;%  -hdac*a
end
figure(3);
xlabel('tag concentration [a]');
ylabel('tagging rate [da/dt]'); % ylabel('tag production rate');
legend(  cellstr(num2str(vs',2)),'Location','Best');
plot([0,50],[0,0],'k-'); hold on;
title('r0 scan, r2=5');
set(gcf,'color','w');