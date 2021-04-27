
% This script produces the violin plots for the condensate model shown in
% Figure 3 and Figure 5. 

% place to save
savePath = 'C:\Shared\Documents\Jordan Looping Model\Revision1\Images\';

% size of the simulation 
tSamples = 100;
nCells = 1000; % we'll track one particular promoter in N isolated stochastic instances, such as in N cells
nR = 21;
ti = 50; %start period

% model parameters 
loopRates = linspace(0,.05,nR); % 
addPol = .01; % can dedimensionalize on this one
losePol = .15;
tSteps = 1e3;
clusterMax =  15; %  max cluster size
% if total PolII is conserved, and bound to diverse sites throughout the
% cell, there is an amount P bound at any given time, distributed across n
% spatially segregated promoters with an average of P/n. If the number of
% total promoters is large, our promoter will rarely have more the average
% available to it.
% A more complete model would have this as a drawn from a random
% distribution of available cluster sizes, which serves only to add spread
% to the on values (as would stochastic transcription models). 

%% start ON
startOn = 1;
promoterPolLog1 = zeros(tSamples,nCells,nR);
for r=1:nR
    promoterPol = startOn*clusterMax*ones(1,nCells);
    e= loopRates(r); % .01;  % looks promising for bistable in e
    tt=0;  
    for t=1:tSteps 
        % promoter gains a Pol molecule
        stoch = rand(nCells,clusterMax+1)  < addPol + e; 
        % its actually faster in Matlab to generate more random numbers
        % than we need (here for Nmax encounters), which saves us having to
        % loop over encouters (even though many loops could exit before
        % Nmax and need fewer total random number draws).  
        for c=1:nCells % would be better to do this without a loop
            stoch(c,promoterPol(c)+2:end) = 0;
        end
        promoterPol(any(stoch,2)) = promoterPol(any(stoch,2))+1;
        promoterPol(promoterPol>clusterMax) = clusterMax;
        % promter loses a pol molecule
        stoch = rand(1,nCells) < losePol;
        promoterPol(stoch) = promoterPol(stoch) - 1;
        promoterPol(promoterPol<0) = 0; 

        if rem(t,tSteps/tSamples)==0
            tt=tt+1;
           promoterPolLog1(tt,:,r) = promoterPol;
        end   
    end
    figure(1); clf; 
    imagesc(promoterPolLog1(:,:,r));
    colorbar; colormap('default');
end

timeAveStON = squeeze(mean(promoterPolLog1(ti:end,:,:),1));
figure(1); clf; imagesc(timeAveStON);
popAveStON = mean(timeAveStON,1);
figure(2); clf; plot(loopRates,popAveStON,'r');

%% start off

startOn = 0;
promoterPolLog0 = zeros(tSamples,nCells,nR);
for r=1:nR
    promoterPol = startOn*clusterMax*ones(1,nCells);
    e= loopRates(r); % .01;  % looks promising for bistable in e
    tt=0; 
    for t=1:tSteps
        stoch = rand(nCells,clusterMax+1)  < addPol + e; % figure(1); clf; imagesc(stoch);
        for c=1:nCells % would be better to do this without a loop
            stoch(c,promoterPol(c)+2:end) = 0;
        end
        promoterPol(any(stoch,2)) = promoterPol(any(stoch,2))+1;
        promoterPol(promoterPol>clusterMax) = clusterMax;

        % promter loses pol
        stoch = rand(1,nCells) < losePol;
        promoterPol(stoch) = promoterPol(stoch) - 1;
        promoterPol(promoterPol<0) = 0; 
        if rem(t,tSteps/tSamples)==0
            tt=tt+1;
           promoterPolLog0(tt,:,r) = promoterPol;
        end   
    end
    figure(1); clf; 
    imagesc(promoterPolLog0(:,:,r));
    colorbar; colormap('default');
end
timeAveStOFF = squeeze(mean(promoterPolLog0(ti:end,:,:),1));
figure(1); clf; imagesc(timeAveStOFF);
popAveStOFF = mean(timeAveStOFF,1);
figure(2); hold on; plot(loopRates,popAveStOFF,'b');

%%
f3 = figure(3); clf;
plot(popAveStON,'r'); hold on;
plot(popAveStOFF,'b'); hold on;
violin(timeAveStON,'bandwidth',.02,'plotMean',false,'faceColor',[1 0 0],'alpha',.5); hold on; 
violin(timeAveStOFF,'bandwidth',.02,'plotMean',false,'faceColor',[0 0 1],'alpha',.5); hold on; 
set(gca,'XTick',1:2:nR,'XTickLabel',num2str(loopRates(1:2:nR)',2));
ylim(round([0,1.05*clusterMax]));
set(gcf,'color','w');

%%


f3 = figure(3); clf;
plot(popAveStON,'r'); hold on;
plot(popAveStOFF,'b'); hold on;
violin(timeAveStON,'bandwidth',.02,'plotMean',false,'faceColor',[1 0 0],'alpha',1); hold on; 
violin(timeAveStOFF,'bandwidth',.02,'plotMean',false,'faceColor',[0 0 1],'alpha',1); hold on; 
set(gca,'XTick',1:2:nR,'XTickLabel',num2str(loopRates(1:2:nR)',2));
ylim(round([0,1.05*clusterMax]));
xlabel('loop rate');
ylabel('Promoter associated PolII');
set(gcf,'color','w');

