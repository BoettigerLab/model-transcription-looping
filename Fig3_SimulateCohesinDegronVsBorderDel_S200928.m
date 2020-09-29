% Stochastic model
% monostable on is still bistable in the stochastic case?
%
% This version accurately follows the determinstic model.

saveFolder = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Images\');
saveName = [saveFolder,'Fig3_Sim.mat'];
%% parameters
k1 = .003; % Rate of conversion from state 0 to state 1. 
k2 = .5; % Rate of conversion from state 1 to state 2. 
tau1 = 20; % Rate of conversion from state 1 to state 0.
tau2 = .05; % Rate of conversion from state 2 to state 1.
hdac = 2; % Rate of tag removal. 
r0 = .01; % Basal tagging rate (state 0)
r1 = .5; % Stimulated tagging rate (state 1) 
r2 =5.0; % Maximally stimulated tagging rate (state 2)
Et= 20; % Maximum tag-transferases associated with promoter

% initial conditions
% start in the active state, explore effect of loss of interaction
%   (symmetric results, the system is hyesteric either way.
Pr_i = 3;
Pr_init = Pr_i*ones(1,Et) ;
ac_init = 40;

% simulation time scales 
t_scale = 1000;
tSteps = 3e6;
nRuns = 3e3;
        
%%
loopRates =[0,6,12];%  [0 5 10];
numLoops = length(loopRates);

%% Run simulation 
if exist(saveName,'file') == 1
    load(saveName,'loopVals','acVals','loopRates');
else
    acVals = cell(numLoops,1);
    loopVals = cell(numLoops,1);
    for r=1:numLoops
        disp(r/numLoops);
        loopRate = loopRates(r);
        acR = [r0 r1 r2]; % acetylation rates for the different enhancer states 
        compTx = false;
        nZymes =Et;
        clear ac_values;
        ds = 100;
        ac_values = cell(nRuns,1); %  zeros(nRuns,tSteps/ds,'uint8');
        loop_events= cell(nRuns,1); % 
        parfor n=1:nRuns
            acR = [r0 r1 r2];
            ac = ac_init; 
            Pr = Pr_init;  % 1, 2, 3
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
    save(saveName,'loopVals','acVals','loopRates');  % loop rates = 0 6 12  
end

%%
av0 =  cat(2,acVals{1}{:}); % enhancer deactivated (equivelantely loop rate=0)
av5 =  cat(2,acVals{2}{:}); % loop rate halved
av10 =  cat(2,acVals{3}{:}); % loop rate high
LoopConds = {av0,av5,av10};
ts = [300 1.4e3 3e4];  % select time 
bDat = cat(1,LoopConds{3}(ts,:),LoopConds{2}(ts,:),LoopConds{1}(ts,:));

figure(4); clf; % as a box plot
hold on; boxplot(bDat');

figure(4); clf; % same data as a violin plot
    violin(bDat','bandwidth',1,'plotMean',false,...
        'faceColor',[.8 .8 .8],'alpha',1,'lineColor',[.7 .7 .7]); 
 hold on;
 plot( mean( bDat,2),'k.','MarkerSize',30);   

figure(3); clf; 
for c=1:3
    subplot(3,1,c);
    violin(LoopConds{4-c}(ts,:)','bandwidth',3,'plotMean',false,...
        'faceColor',[.8 .8 .8],'alpha',1,'lineColor',[.7 .7 .7]); 
    hold on;
    med = mean( double(LoopConds{4-c}(ts,:)),2 );
    plot(1:3,med,'k.','MarkerSize',30);
end
