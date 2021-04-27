
% condesnate model parameter exploration
% This script scans over the parameter space of the condensate model
% This produces the heat maps shown in Figure 3. 

% size of the simulation
tSamples = 100; % compressed for speed
tSteps = 1e3;
nCells = 500; %
nR = 51;
ti = 50;

% parameter scans
nP1 = 10; % number of parameters to test
nP2 = 20; % number of parameters to test
loopRates = linspace(0,.2,nR)'; % 
clusterMaxes = ceil(logspace(log10(1),log10(30),nP1)); %  1:2:2*nP1;
offRates = logspace(log10(.01),log10(1),nP2); %  linspace(0.01,1,nP2);
startOn = 0;

yResults = zeros(nP1,nP2,nR);
gainMatrixH = zeros(nP1,nP2); % steepness in Hill form
switchMatrixH = zeros(nP1,nP2); % km hill
maxMatrixH = zeros(nP1,nP2); % 
parfor p1=1:nP1
    for p2 = 1:nP2
        % model parameters      
        addPol = .01; % can dedimensionalize on this one
        losePol = offRates(p2); %#ok<*PFBNS> % .15;
        clusterMax =  clusterMaxes(p1); %  5;% 
        promoterPolLog = zeros(tSamples,nCells,nR);
        for r=1:nR
            promoterPol = startOn*clusterMax*ones(1,nCells);
            e= loopRates(r); 
            % start dynamic simulation
            tt=0;
            for t=1:tSteps 
                % Gain event
                stoch = rand(nCells,clusterMax+1)  < addPol + e; % 10% time
                for c=1:nCells   % would be great to lose the loop here. 
                    stoch(c,promoterPol(c)+2:end) = 0;  % 85 time
                end
                promoterPol(any(stoch,2)) = promoterPol(any(stoch,2))+1;
                promoterPol(promoterPol>clusterMax) = clusterMax;
                % Loss event
                stoch = rand(1,nCells) < losePol;
                promoterPol(stoch) = promoterPol(stoch) - 1;
                promoterPol(promoterPol<0) = 0; 
                % periodically system state in log
                if rem(t,tSteps/tSamples)==0
                    tt=tt+1;
                   promoterPolLog(tt,:,r) = promoterPol;
                end   
            end
        end
        % compute curve
         timeAve = squeeze(mean(promoterPolLog(ti:end,:,:),1));
         popAve = mean(timeAve,1); 
        % fit popAve with Hill
        x = loopRates;
        y = popAve'; % convert to a column vector
        yResults(p1,p2,:) = y;        
        hillEqn = 'a*x^c/(b^c+x^c)';
        fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[1.2*clusterMax,max(x),inf],...
                   'StartPoint',[clusterMax,mean(x),4]);
        fitResult = fit(x,y,hillEqn,fo);
        gainMatrixH(p1,p2) = fitResult.c;
        switchMatrixH(p1,p2) = fitResult.b;
        maxMatrixH(p1,p2) = fitResult.a;
    end
end

%% in Hill form
f2=figure(2); clf;
subplot(1,3,1); imagesc(log2(gainMatrixH(:,1:end-1)));
set(gca,'Ytick',1:nP1,'YtickLabel',clusterMaxes);
set(gca,'Xtick',1:5:nP2,'XtickLabel',offRates(1:5:end)); colorbar;
title('gain matrix'); colorbar; xlabel('off rate'); ylabel('cMax');
caxis([log2(0.1),log2(10)]);
subplot(1,3,2); imagesc(log2(switchMatrixH(:,1:end-1))); 
set(gca,'Ytick',1:nP1,'YtickLabel',clusterMaxes);
set(gca,'Xtick',1:5:nP2,'XtickLabel',offRates(1:5:end)); colorbar;
title('switch matrix'); colorbar; xlabel('off rate'); ylabel('cMax'); 
caxis(log2([min(loopRates),max(loopRates)]));
subplot(1,3,3); imagesc(log2(maxMatrixH(:,1:end-1)));
set(gca,'Ytick',1:nP1,'YtickLabel',clusterMaxes);
set(gca,'Xtick',1:5:nP2,'XtickLabel',offRates(1:5:end)); colorbar;
title('max matrix'); colorbar;
xlabel('off rate'); ylabel('cMax');
caxis([log2(1),log2(30)]);
set(gcf,'color','w'); 

figure(1); clf; plot(loopRates,squeeze(yResults(10,12,:))); hold on;
plot(loopRates,squeeze(yResults(5,end,:))); hold on;
plot(loopRates,squeeze(yResults(3,end,:))); hold on;
 
%% Save Result
% update file path
% note, SetFigureSavePath and SaveFigure require Boettiger-lab
% matlab-functions.  If these functions are not installed you may simply
% save the figures using matlab's built in save options. 
saveFolder = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Revision1\Images\'); 
SaveFigure(f2,'name','condensateParScan_log','formats',{'eps','png','fig'},'overwrite',false);




    
    