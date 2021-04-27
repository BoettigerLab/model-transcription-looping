%% simulation of two different promoters under same 2x decrease in contact
% achieved with two different values of losePol (factor of 2 diff)
% this is the condensate version of the model

savePath = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Revision1\Images\');

% size of the simulation 
tSamples = 100; % total number of samples to grab in time 
nCells = 1000; % we'll track one particular promoter in N isolated stochastic instances, such as in N cells
ti = 50; %start period - we'll wait this long for things to equilibrate

% model parameters 
addPol = .01; % can dedimensionalize on this one
tSteps = 1e3;
clusterMax =  15; %  max cluster size
startOn = 0;
e = .05;
losePol_array = [.15,.30];
loopRates = [.015,.030];
nR = 4; % total simulation conditions
promoterPolLog0 = zeros(tSamples,nCells,nR);
r = 0;
for a=1:length(losePol_array) % loop over the 2 promoters
    for b=1:length(loopRates) % loop over the 2 different loop rates
        r=r+1;
        % update variable parameters
        promoterPol = startOn*clusterMax*ones(1,nCells);
        e =  loopRates(b); % .01;  % looks promising for bistable in e
        losePol = losePol_array(a); 
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
end
%% plotting
timeAveStOFF = squeeze(mean(promoterPolLog0(ti:end,:,:),1));
popAveStOFF = mean(timeAveStOFF,1);
o = [1,3,2,4];
names = {'P_1 1x','P_1 2x','P_2 1x','P_2 2x'};
f3 = figure(3); clf;
violin(timeAveStOFF(:,o),'bandwidth',.02,'plotMean',false,'faceColor',[0 0 1],'alpha',1); hold on;  % can't save as eps with alpha
ylim(round([0,1.05*clusterMax]));
set(gcf,'color','w');
ylabel('PolII at Promoter')
hold on; plot(1:length(names),nanmedian(timeAveStOFF(:,o)),'r.','MarkerSize',10);
hold on; plot(1:length(names),nanmean(timeAveStOFF(:,o)),'m+','MarkerSize',4);
set(gca,'xTick',1:length(names),'XTickLabels',names(o));
