%% simulations at two different values of losePol (factor of 2 diff)


savePath = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Revision1\Images\');

% size of the simulation 
tSamples = 100; % total time samples to save
nCells = 1000; % we'll track one particular promoter in N isolated stochastic instances, such as in N cells
ti = 50; %start period

% model parameters 
addPol = .01; % can dedimensionalize on this one
tSteps = 1.5e3;
clusterMax =  15; %  max cluster size
startOn = 1;
e = .05;
losePol_array = [.22];
loopRates = [0,.015,.030];
nR = length(losePol_array)  * length(loopRates); % total simulations
promoterPolLog0 = zeros(tSamples,nCells,nR);
%% simulation
r = 0;
for a=1:length(losePol_array) % loop over promoter conditions (fixed/common)
    for b=1:length(loopRates) % loop over different loop rates
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
timeAveStOFF = squeeze(mean(promoterPolLog0(ti:end,:,:),1));
popAveStOFF = mean(timeAveStOFF,1);

%% plotting

f3 = figure(3); clf;
violin(timeAveStOFF,'bandwidth',.02,'plotMean',false,'faceColor',[0 0 1],'alpha',.5); hold on; 
ylim(round([0,1.05*clusterMax]));
set(gcf,'color','w');

timeAveEarly = squeeze(mean(promoterPolLog0(9:11,:,:),1));
timeAveLate = squeeze(mean(promoterPolLog0(80:end,:,:),1));
tempComp = [timeAveEarly,timeAveLate];
f4 = figure(4); clf;
violin(tempComp,'bandwidth',.02,'plotMean',false,'faceColor',[0 0 1],'alpha',1); hold on; 
ylim(round([0,1.05*clusterMax]));
set(gcf,'color','w');
names = {'Early Del En','Early 1/2x','Early 1x','Late Del En','Late 1/2x','Late 1x'};
set(gca,'xTick',1:length(names),'XTickLabels',names);
hold on; plot(1:length(names),nanmedian(tempComp),'r.','MarkerSize',10);
hold on; plot(1:length(names),nanmean(tempComp),'m+','MarkerSize',4);



