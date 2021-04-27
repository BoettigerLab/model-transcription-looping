%% Deterministic Model Analysis
% PTM model derivation of analytic results in the ODE continuum limit. 

% states
% a0 <=> a1 <=> a2  ... <=> an  % accumulating amount of tag 'a'
% Pr <=> Pr-E1 <=> Pr-E2
SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Images\');
syms E0 E1 E2 Et tau k1 k2 r0 r1 r2 hdac a;
consE = Et -(E0 + E1 + E2);
dE0 = -k1*a*E0 + tau*E1;
dE1 =  k1*a*E0 + tau*E2 -k2*a*E1 -tau*E1;
dE2 =  k2*a*E1 - tau*E2;
nullE = solve([dE0==0,dE1==0,dE2==0,consE==0],[E0,E1,E2]);
nullE.E0
nullE.E1
nullE.E2
da = r0*E0 + r1*E1 + r2*E2 - hdac*a;

simplify( r0*nullE.E0 + r1*nullE.E1 + r2*nullE.E2)

%%

% high a highly separated parameters
 %  % syms a r0 r1 r2 k1 k2 tau hdac loopRate
% note, since the stochastic model moves in integer units of ac, the
% absolute value of ac matters, its not just a free scaling parameter. 
k1 = .01; % higher steepens the curve, pulling the Km to lower values. 
k2 = .4; % higher steepens the curve, pulling the Km to lower values
%   also k2/k1 large -> more sigmoidal k1>k2 not sigmoidal
tau = 1; % just scales the ks

hdac = 17; % just scales the rs. 
r0 = 1; % .7; % 3/1000/.005; % sets behavior at ac=0, effects loopRate contribution
r1 = 5; % not super sensitive 
r2 =50; % large values for low ks pushes Km to higher values
loopRate = 28; % .65 ;% .95; % the linear shift to the sigmoid. % 20% change, 3.5 to 4.5  
%da/dt = (r0+r1*k1/tau*a + r2*(k1*k2/tau^2)*a^2)/(1 + k1/tau*a + k1*k2/tau^2*a^2) - hdac*a+loopRate;
Et= 10;
a = 0:0.1:60; 
p = (Et*(k1*k2*r2*a.^2 + k1*r1*a*tau + r0*tau^2))./(k1*k2*a.^2 + k1*a*tau + tau^2) +loopRate;
figure(1); clf; plot(a,p); hold on; plot(a,hdac*a); pause(.1);
%%
syms a
nSteps = 100;
loopRates = linspace(0,45,nSteps);
lows = zeros(nSteps,1);
highs = zeros(nSteps,1);
for i=1:nSteps
    dadt = (Et*(k1*k2*r2*a.^2 + k1*r1*a*tau + r0*tau^2))./(k1*k2*a.^2 + k1*a*tau + tau^2) +loopRates(i) -a*hdac;
    soln = vpasolve(dadt==0,a,[0,60]);
    lows(i) = min(soln);
    highs(i) = max(soln);
end

f2 = figure(2); clf;
plot(loopRates,lows,'bo'); hold on;
plot(loopRates,highs,'r.');
plot(loopRates,lows,'b-'); hold on;
plot(loopRates,highs,'r-');
xlabel('contact frequency');
ylabel('transcription rate');
set(gcf,'color','w');

SaveFigure(f2,'name','ex_hystersis','formats',{'png','eps'},'overwrite',true);



%%


% high ac highly separated parameters
 %  % syms a r0 r1 r2 k1 k2 tau hdac loopRate
% note, since the stochastic model moves in integer units of ac, the
% absolute value of ac matters, its not just a free scaling parameter. 
k1 = .005; % higher steepens the curve, pulling the Km to lower values. 
k2 = .5; % higher steepens the curve, pulling the Km to lower values
%   also k2/k1 large -> more sigmoidal k1>k2 not sigmoidal
tau = 1; % just scales the ks

hdac = 13; % just scales the rs. 
r0 = 1; % .7; % 3/1000/.005; % sets behavior at ac=0, effects loopRate contribution
r1 = 5; % not super sensitive 
r2 =50; % large values for low ks pushes Km to higher values
loopRate = 25; % .65 ;% .95; % the linear shift to the sigmoid. % 20% change, 3.5 to 4.5  
%da/dt = (r0+r1*k1/tau*a + r2*(k1*k2/tau^2)*a^2)/(1 + k1/tau*a + k1*k2/tau^2*a^2) - hdac*a+loopRate;
Et= 10;
a = 0:0.1:60; 
p = (Et*(k1*k2*r2*a.^2 + k1*r1*a*tau + r0*tau^2))./(k1*k2*a.^2 + k1*a*tau + tau^2) +loopRate;
figure(1); clf; plot(a,p); hold on; plot(a,hdac*a); pause(.1);
%%
syms a
nSteps = 100;
loopRates = linspace(0,45,nSteps);
lows = zeros(nSteps,1);
highs = zeros(nSteps,1);
for i=1:nSteps
    dadt = (Et*(k1*k2*r2*a.^2 + k1*r1*a*tau + r0*tau^2))./(k1*k2*a.^2 + k1*a*tau + tau^2) +loopRates(i) -a*hdac;
    soln = vpasolve(dadt==0,a,[0,60]);
    lows(i) = min(soln);
    highs(i) = max(soln);
end

f2 = figure(2); clf;
plot(loopRates,lows,'bo'); hold on;
plot(loopRates,highs,'r.');
plot(loopRates,lows,'b-'); hold on;
plot(loopRates,highs,'r-');
xlabel('contact frequency');
ylabel('transcription rate');
set(gcf,'color','w');
%% 
SaveFigure(f2,'name','ex_hystersis_memory','formats',{'png','eps'},'overwrite',true);


%% phase portraits

%da/dt = (r0+r1*k1/tau*a + r2*(k1*k2/tau^2)*a^2)/(1 + k1/tau*a + k1*k2/tau^2*a^2) - hdac*a+loopRate;

eqn = 'da/dt = E_T*(r_2*a.^2 + r_1*a*(\tau/k_2) + r_0*\tau^2/(k_1*k_2)))./(a.^2 + a*(\tau/k_2) + \tau^2/(k_1*k_2)) +loopRate - d_{hdac}*a';
figure(1); clf; title(eqn); set(gcf,'color','w');

%
% slope at a = Km = tau/sqrt(k_1*k_2)

%%

% high ac highly separated parameters
 %  % syms a r0 r1 r2 k1 k2 tau hdac loopRate
% note, since the stochastic model moves in integer units of ac, the
% absolute value of ac matters, its not just a free scaling parameter. 
k1 = .5; % higher steepens the curve, pulling the Km to lower values. 
k2 = .2; % higher steepens the curve, pulling the Km to lower values
%   also k2/k1 large -> more sigmoidal k1>k2 not sigmoidal
tau = 10; % just scales the ks

hdac = 1.8; % just scales the rs. 
r0 = .00; % .7; % 3/1000/.005; % sets behavior at ac=0, effects loopRate contribution
r1 = 10; % not super sensitive 
r2 =50; % large values for low ks pushes Km to higher values
loopRate = 0; % .65 ;% .95; % the linear shift to the sigmoid. % 20% change, 3.5 to 4.5  
%da/dt = (r0+r1*k1/tau*a + r2*(k1*k2/tau^2)*a^2)/(1 + k1/tau*a + k1*k2/tau^2*a^2) - hdac*a+loopRate;
Et= 3;
a = 0:0.1:180; 
p = (Et*(k1*k2*r2*a.^2 + k1*r1*a*tau + r0*tau^2))./(k1*k2*a.^2 + k1*a*tau + tau^2) +loopRate;
% limiting max Et (r0+r1+r2) bounds. r2 behavior dominates at large a. r0 behavior dominates at small a.
%     so r0 < r2  in order for bistability (to cross the linear a*hdac expression it must be positive sigmoid.  
%     can still be bistable for r1 > r2, but not r1 >> r2 becomes monostable  
figure(1); clf; plot(a,p); hold on; plot(a,hdac*a); pause(.1);

syms a
dadt = Et*(r0+r1*k1/tau*a + r2*(k1*k2/tau^2)*a^2)/(1 + k1/tau*a + k1*k2/tau^2*a^2) +loopRate - hdac*a;

% dp/da = hdac; % parameters where derivative of production with respect to
% a equals hdac at the position where d^2p/dt^2 = 0.
% AND where the linear regime is as large as possible 
%    Taylor expand the function around the point at which d^2p/dt^2 = 0.

vpasolve(dadt == 0,a)
tau/sqrt(k1*k2)
%%
syms E0 E1 E2 Et tau k1 k2 r0 r1 r2 hdac a;
syms a;
p = (r2*a.^2 + r1*a*tau/k2 + r0*tau^2/(k1*k2))./(a.^2 + a*tau/k2 + tau^2/(k1*k2));
dp = diff(p,a);
subs(dp,a,tau/sqrt(k1*k2))


syms E0 E1 E2 Et tau k1 k2 r0 r1 r2 hdac a;
syms a;
p = (a.^2 )./(a.^2 + a*tau/k2 + tau^2/(k1*k2));
dp = diff(p,a);
subs(dp,a,tau/sqrt(k1*k2))


(2*Et*(k1*tau + 2*a*k1*k2)^2*(k1*k2*r2*a^2 + k1*r1*a*tau + r0*tau^2))/(k1*k2*a^2 + k1*a*tau + tau^2)^3 ...
+ (2*Et*k1*k2*r2)/(k1*k2*a^2 + k1*a*tau + tau^2) ...
- (2*Et*(k1*r1*tau + 2*a*k1*k2*r2)*(k1*tau + 2*a*k1*k2))/(k1*k2*a^2 + k1*a*tau + tau^2)^2 ...
- (2*Et*k1*k2*(k1*k2*r2*a^2 + k1*r1*a*tau + r0*tau^2))/(k1*k2*a^2 + k1*a*tau + tau^2)^2
 

vpasolve(ddp==0,a)
