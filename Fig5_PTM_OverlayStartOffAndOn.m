
% load data from start high and start low simulations
loopDataFolder = 'U:\Data\Analysis\2020-04-23_JordanLoopSims\'; % update for your save path 
load([loopDataFolder,'loopRatesVsExpStoch.mat'],'acVals2'); % start low
load([loopDataFolder,'loopRatesVsExpStoch_startON.mat'],'acVals3'); % start high

numLoops = length(acVals2);
[nRuns, tSteps] = size(acVals2{1});
lowDist = zeros(numLoops,nRuns);
highDist= zeros(numLoops,nRuns);
for r=1:numLoops
    endValsLow = acVals2{r}(:,end);
    endValsHigh = acVals3{r}(:,end);
    lowDist(r,:) = endValsLow;
    highDist(r,:) = endValsHigh;
end

%% high
f3 = figure(3); clf;
subsample = 1:5:numLoops-10;
% % opitionally overlay violin plots of start low
% violin(lowDist(subsample,:)','bandwidth',.5,'plotMean',false,'faceColor',[.5 .5 1],'alpha',1); hold on;
% violin plots of start high
violin(highDist(subsample,:)','bandwidth',.5,'plotMean',false,'faceColor',[1 .5 .5],'alpha',1); hold on; 

% overlay medians
lowMed = median(lowDist(subsample,:),2);
highMed = median(highDist(subsample,:),2);
plot(lowMed,'b'); hold on;
plot(highMed,'r');
ylim([0,70]);
set(gca,'xTickLabel',loopRates(subsample));
set(gcf,'color','w');
xlabel('loop rate');
ylabel('transcription rate');

%% Temporal topographic plot
f10 = figure(10); clf;
medLow = zeros(numLoops,1);
medHigh = zeros(numLoops,1);
ts=10:200:3000;
blueMap = GetColorMap('blue',length(ts));
redMap = GetColorMap('red',length(ts));
for i=1:length(ts)
    t=ts(i);
    for r=1:numLoops
       medLow(r) = median(acVals2{r}(:,t));
        medHigh(r) = median(acVals3{r}(:,t));
    end
    plot(loopRates,medLow,'-','linewidth',1,'color',blueMap(length(ts)-i+1,:));  hold on;
    plot(loopRates,medHigh,'-','linewidth',1,'color',redMap(length(ts)-i+1,:)); 
    set(gcf,'color','w');
end
xlabel('contact frequency');
ylabel('median transcription');
set(gca,'color',[.8 .8 .8]);
set(gcf,'color','w');
