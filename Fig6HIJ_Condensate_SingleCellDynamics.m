
savePath = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Revision1\Images\');
%% single cell correlation of contact and transcription bursts

% size of the simulation 
tSamples = 1e3; % total number of samples 
nCells = 1000; % we'll track one particular promoter in N isolated stochastic instances, such as in N cells
ti = 50; % start period -- wait this long after initializing

% model parameters 
addPol = .01; % can dedimensionalize on this one
% losePol = .15;
tSteps = 1e3;
clusterMax =  15; %  max cluster size
startOn = 0;
losePol_array = .15; % [.15,.30];
loopRates = linspace(0,.05,21); % 
nR = length(loopRates);
promoterPolLog0 = zeros(tSamples,nCells,nR);
loopEvents = zeros(tSamples,nCells,nR);
r = 0;
for a=1:length(losePol_array)
    for b=1:length(loopRates)
        r=r+1;
        % update variable parameters
        promoterPol = startOn*clusterMax*ones(1,nCells);
        e =  loopRates(b); % .01;  % looks promising for bistable in e
        losePol = losePol_array(a);
        tt=0;
        for t=1:tSteps
            r_stoch = rand(nCells,clusterMax+1) ;
            stoch = r_stoch< addPol + e; % figure(1); clf; imagesc(stoch);
            didLoop = r_stoch(:,1) < e; % the raw looping rate (don't want to call promoters with condensate as more freq loopers, just more frequent binders)
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
               loopEvents(tt,:,r) = didLoop;
            end   
        end
        figure(1); clf; 
        imagesc(promoterPolLog0(:,:,r));
        colorbar; colormap('default');
    end
end
timeAveStOFF = squeeze(mean(promoterPolLog0(ti:end,:,:),1));
popAveStOFF = mean(timeAveStOFF,1);

%% 
f3 = figure(3); clf;
violin(timeAveStOFF,'bandwidth',.02,'plotMean',false,'faceColor',[0 0 1],'alpha',.5); hold on; 
ylim(round([0,1.05*clusterMax]));
set(gcf,'color','w');
ylabel('PolII at Promoter')
%%
r=11;
y1= loopEvents(:,c,r); % time x cell x looprate
y2 = promoterPolLog0(:,c,r);
figure(1); clf; subplot(1,2,1); plot(y1);
subplot(1,2,2); plot(y2);
%%
y1= loopEvents(:,:,r); % time x cell x looprate
y2 = promoterPolLog0(:,:,r);
binsWithContact = y1(:)>.5;
binsWithTx = y2(:) > clusterMax-1;
[odds, cI] = OddsRatioCI(binsWithContact,binsWithTx,'iters',10)

%% larger bin sizes
bn  = 2;
odds = zeros(nR,1);
cI = zeros(nR,2);
for r=1:nR
    binsWithContact = imresize( uint8(255*loopEvents(:,:,r)),[tSamples/bn,nCells]);
    binsWithTx = imresize( uint16(255*promoterPolLog0(:,:,r)),[tSamples/bn,nCells]);
    [odds(r), cI(r,:)] = OddsRatioCI(binsWithContact(:)>1,binsWithTx(:)>clusterMax*20,'iters',10);
end

figure(1); clf;
r = 8:nR;
ploterr(loopRates(r),odds(r),[],{cI(r,1),cI(r,2)},'.','color','k');
hold on; plot([.012,.052],[1,1],'k--');
xlim([.014,.053 ]);
ylim([0,2]);
xlabel('loop rate (arbitary units)');
ylabel('Odds Ratio');
%%

%% Fig. 5E time from burst
% loops = double(cat(2,loopVals{4/.2}{:})) > 0;
% save([savePath,'fig5condensate_sim.mat'],'loopRates','loopEvents','promoterPolLog0'); % for exact reproduction since its a stochastic sim;
r =12;
bn = 1;
binsWithContact = imresize( uint8(255*loopEvents(:,:,r)),[tSamples/bn,nCells]) > 1;
binsWithTx = imresize( uint16(255*promoterPolLog0(:,:,r)),[tSamples/bn,nCells]) > clusterMax*20;
  
w=100; % 40
nC_low =  400:428; % 350:450; % %% arbitrary subset of cells
nC_high = 1:1000; 
nCs =  nC_high; % nC_low; % 
nC = length(nCs);
t1 =250:300; % (tSteps/bn)-2*w;
burstCentered = cell(nC*length(t1),1);
k=0;
for c=nCs
    for t=t1
        if binsWithTx(t,c) % has expression burst
            k=k+1;
            burstCentered{k} = binsWithContact(t-w+1:t+w,c)';
        end
    end
end
burstCentered = cat(1,burstCentered{:});
m =  nanmean( burstCentered,1);
s1 =  m + nanstd( burstCentered,1)/sqrt(size(burstCentered,1));
s2 =  m - nanstd( burstCentered,1)/sqrt(size(burstCentered,1));
f1 = figure(1); clf; plot(-w+1:w,1-m ); hold on;
plot(-w+1:w,1-s1,'--');
plot(-w+1:w,1-s2,'--');
title(['N bursts = ',num2str(k)]);
ylabel('Ave. Distance 0=contact 1=no contact')
set(gcf,'color','w');
xlabel('time relative to burst');
ylim([0.94,1.01]);
