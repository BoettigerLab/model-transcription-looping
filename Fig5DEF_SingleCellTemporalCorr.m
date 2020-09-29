% Stochastic model
% monostable on is still bistable in the stochastic case?
%
% This version accurately follows the determinstic model.

saveFolder = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Images\');
saveName = [saveFolder,'loopResultsStartLow.mat'];
%% Parameters
k1 = .003; % Rate of conversion from state 0 to state 1. 
k2 = .5; % Rate of conversion from state 1 to state 2. 
tau1 = 20; % Rate of conversion from state 1 to state 0.
tau2 = .05; % Rate of conversion from state 2 to state 1.
hdac = 2; % Rate of tag removal. 
r0 = .01; % Basal tagging rate (state 0)
r1 = .5; % Stimulated tagging rate (state 1) 
r2 =5.0; % Maximally stimulated tagging rate (state 2)
Et= 20; % Maximum tag-transferases associated with promoter

% initial conditions: start off, unstimulated enzyme, no promoter tags
ac_init = 0; 
Pr_init = 1;


% simulation time scales  
t_scale = 1000;
tSteps = 3e5;
nRuns = 3e3;
    %%
loopRates = 0:.2:20;
numLoops = length(loopRates);

if exist(saveName,'file')==0
acVals = cell(numLoops,1);
loopVals = cell(numLoops,1);
for r=1:numLoops
    disp(r/numLoops);
    loopRate = loopRates(r);
    acR = [r0 r1 r2];  % acetylation rates for the different enhancer states 
    compTx = false;
    nZymes =Et;
    clear ac_values;
    ds = 100; % downsample to save storage spac
    ac_values = cell(nRuns,1); % 
    loop_events= cell(nRuns,1); % 
    parfor n=1:nRuns
        acR = [r0 r1 r2];
        ac = ac_init; % 0 ... inf
        Pr = Pr_init*ones(1,nZymes);  % 1, 2, 3
        ac_values{n} = zeros(tSteps/ds,1,'uint8');
        loop_events{n} = zeros(tSteps/ds,1,'uint8');
        loopCount = 0;
        for t=1:tSteps % straight discreet time simulation
            % enhancer mediated actyl
            if rand*t_scale < loopRate
               ac = ac + 1; 
               loopCount = loopCount+1;
            end
            % lose actyl
            if (rand*t_scale < hdac*ac) && ac > 1
                ac = ac-1;
            end
            % handle promoter state
            for e=1:nZymes
                % gain actyl
                if rand*t_scale < acR(Pr(e)) % rate depends on enzyme state
                    ac = ac + 1;
                end
                if Pr(e)==1
                    if rand*t_scale < ac*k1
                        Pr(e)=2;
                    end
                elseif Pr(e)==2 
                    if rand*t_scale < ac*k2
                        Pr(e) = Pr(e)+1;
                    end
                    if rand*t_scale < tau1
                        Pr(e) = Pr(e)-1;
                    end
                elseif Pr(e) == 3
                    if rand*t_scale < tau2
                        Pr(e) = Pr(e) - 1;
                    end
                end
            end
            if rem(t,100) == 0 % downample to save memory
                tt=t/100;
                ac_values{n}(tt) = ac;
                loop_events{n}(tt) = loopCount; % record number of loops in last sampling interval
                loopCount = 0; % reset to zero
            end
        end
    end
    acVals{r} = ac_values; % 
    loopVals{r} = loop_events;
end
    save(saveName,'loopVals','acVals','loopRates');
else
    load(saveName,'loopVals','acVals','loopRates');
end
%%

    
%% Supp Fig. 5B: percent time bursting vs time in contact
% 500 cells, bar graphs
hiCell_loops = double(cat(2,loopVals{4/.2}{:}));
hiCell_expr = double(cat(2,acVals{4/.2}{:}));
hiCell_loops = hiCell_loops(1501:2000,100:900);
hiCell_expr = hiCell_expr(1501:2000,100:900);
hasBurst = hiCell_expr > 20;
m = nan(5,1);
s=nan(5,1);
for i=1:5
    lr = 5-i;
    m(i) = mean( hasBurst(hiCell_loops==lr))*100;
    s(i) = std( hasBurst(hiCell_loops==lr)) / sqrt( length(hasBurst(hiCell_loops==lr)) )*100;
end
f2 = figure(2); clf; 
bar(m); hold on;
ploterr(1:5,m,[],s,'.'); 
ylim([0,5]); ylabel('percent time bursting');
xlabel('loop rate bin (high to low)  ~ave distance bin');

%% Supp. Fig 5C. high sampled
x=round(4.8/.2);
hiCell_loops = double(cat(2,loopVals{x}{:}));
hiCell_expr = double(cat(2,acVals{x}{:}));
hasBurst = hiCell_expr > 20;
m = nan(5,1);
s=nan(5,1);
for i=1:5
    lr = 5-i;
    m(i) = mean( hasBurst(hiCell_loops==lr))*100;
    s(i) = std( hasBurst(hiCell_loops==lr)) / sqrt( length(hasBurst(hiCell_loops==lr)) )*100;
end
f2 = figure(2); clf; 
bar(m); hold on;
ploterr(1:5,m,[],s,'.'); 
ylim([0,5]); ylabel('percent time bursting');
xlabel('loop rate bin (high to low)  ~ave distance bin');

%% Supp Fig5 E boxplot comparison
hiCell_loops = double(cat(2,loopVals{4/.2}{:}));
hiCell_expr = double(cat(2,acVals{4/.2}{:}));
hasBurst = hiCell_expr > 20;
hasLoop = hiCell_loops>0;
f1 = figure(1); clf;
on = double(hasLoop); on(~hasBurst) = NaN;
off= double(hasLoop); off(hasBurst) = NaN;
% A small subset of the data: Supp Fig 5E.
cellSubset1 = 1:180;   timeSubset1 = 800:1230;  cellSubset2 = 1:180 ;    timeSubset2 = 800:1230; % randperm(size(on,1),100);
onDist = 1- nanmean(on(timeSubset1,cellSubset1),1);
offDist = 1-nanmean(off(timeSubset2,cellSubset2),1);
BoxPlot2D({onDist,offDist},[],'wisker',.9); hold on;  
[p,h] = ranksum(onDist,offDist);
disp(p);
ylim([0,1]); xlim([0,3]); xlabel('ON frames     OFF frames');
ylabel('average distance 0=contact 1=no contact');
   
%% Supp Fig5  F:  boxplot comparison
hiCell_loops = double(cat(2,loopVals{4/.2}{:}));
hiCell_expr = double(cat(2,acVals{4/.2}{:}));
hasBurst = hiCell_expr > 20;
hasLoop = hiCell_loops>0;
f1 = figure(1); clf;
on = double(hasLoop); on(~hasBurst) = NaN;
off= double(hasLoop); off(hasBurst) = NaN;
% A larger subset of the data: Supp Fig 5F.
cellSubset1 =1:3000; timeSubset1 =1:3000; cellSubset2 =1:3000; timeSubset2 =1:3000; 
onDist = 1- nanmean(on(timeSubset1,cellSubset1),1);
offDist = 1-nanmean(off(timeSubset2,cellSubset2),1);
BoxPlot2D({onDist,offDist},[],'wisker',.9); hold on;  
[p,h] = ranksum(onDist,offDist);
disp(p);
ylim([0,1]); xlim([0,3]); xlabel('ON frames     OFF frames');
ylabel('average distance 0=contact 1=no contact');
  


%% Fig. 5E time from burst
loops = double(cat(2,loopVals{4/.2}{:})) > 0;
expr = double(cat(2,acVals{4/.2}{:}));    
w=40;
nC = 100; %
t1 =1500:1800;
t2 = 1000:2500;
loopCentered = cell(nC*300,1);
k=0;
for c=1:nC
    for t=t1
        if expr(t,c) > 20
            k=k+1;
            loopCentered{k} = loops(t-w+1:t+w,c)';
        end
    end
end
loopCentered = cat(1,loopCentered{:});
m =  nanmean( loopCentered,1);
s1 =  m + nanstd( loopCentered,1)/sqrt(size(loopCentered,1));
s2 =  m - nanstd( loopCentered,1)/sqrt(size(loopCentered,1));
f1 = figure(1); clf; plot(-w+1:w,1-m ); hold on;
plot(-w+1:w,1-s1,'--');
plot(-w+1:w,1-s2,'--');
title(['N bursts = ',num2str(k)]);
ylabel('Ave. Distance 0=contact 1=no contact')
set(gcf,'color','w');
xlabel('time relative to burst');
ylim([0.6,.8]);

%% Fig. 5F time from burst
loops = double(cat(2,loopVals{4/.2}{:})) > 0;
expr = double(cat(2,acVals{4/.2}{:}));    
w=40;
nC =  3000 % 100; %
t1 =1500:1800;
t2 = 1000:2500;
loopCentered = cell(nC*300,1);
k=0;
for c=1:nC
    for t=t1
        if expr(t,c) > 20
            k=k+1;
            loopCentered{k} = loops(t-w+1:t+w,c)';
        end
    end
end
loopCentered = cat(1,loopCentered{:});
m =  nanmean( loopCentered,1);
s1 =  m + nanstd( loopCentered,1)/sqrt(size(loopCentered,1));
s2 =  m - nanstd( loopCentered,1)/sqrt(size(loopCentered,1));
f1 = figure(1); clf; plot(-w+1:w,1-m ); hold on;
plot(-w+1:w,1-s1,'--');
plot(-w+1:w,1-s2,'--');
title(['N bursts = ',num2str(k)]);
ylabel('Ave. Distance 0=contact 1=no contact')
set(gcf,'color','w');
xlabel('time relative to burst');
ylim([0.6,.8]);

    


