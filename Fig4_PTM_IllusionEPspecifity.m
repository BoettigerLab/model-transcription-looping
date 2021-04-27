% Stochastic model
% monostable on is still bistable in the stochastic case?
%
% This version accurately follows the determinstic model.

saveFolder = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Images\');

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
theta = 0; 
Pr_i = 1;

% two different intrinsic promoter hdac affinities
hdacs = [2,4]; 
% loop rates, before and after perturbation, for each promoter
loopRates = [5 10]; 

% number of simulations, length of simulation, time scale of sampling
tSteps = 3e5; % length of simulation
nRuns = 3e3; % number of replicates
t_scale = 1000; % time scale of sampling

saveName = [saveFolder,'loopResultsStartLowEPSpec.mat']; 

%% Run simulation 
if exist(saveName,'file') == 1
    load(saveName,'loopVals','acVals','loopRates');
else

%% Run simulation
    numLoops = length(loopRates);
    acVals = cell(numLoops,2);
    loopVals = cell(numLoops,2);
    for h=1:2
        hdac = hdacs(h); 
        for r=1:numLoops
            disp(r/numLoops);
            loopRate = loopRates(r);
            % this version generally falls off and stays off
            ac_init = 2*theta; % initial
            Pr_init = Pr_i;
            acR = [r0 r1 r2]; % [r0,r1,r2] % acetylation rates for the different enhancer states 
            compTx = false;
            nZymes =Et;
            clear ac_values;
            ds = 100;
            ac_values = cell(nRuns,1); %  zeros(nRuns,tSteps/ds,'uint8');
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
            acVals{r,h} = ac_values; % 
            loopVals{r,h} = loop_events;
        end
    end
    save(saveName,'loopVals','acVals','loopRates');
end

%%
p1_lo =  cat(2,acVals{1,1}{:}); % low hdac, low loop
p1_hi =  cat(2,acVals{2,1}{:}); % low hdac, hi loop
p2_lo =  cat(2,acVals{1,2}{:}); % higher hdac, low loop
p2_hi =  cat(2,acVals{2,2}{:}); % higer hdac, hi loop
t = 2000; x = 0:2:70;
vdists = [p1_lo(t,:); p2_lo(t,:);  p1_hi(t,:); p2_hi(t,:)]';

f4=figure(4); clf;
violin(vdists,'bandwidth',1,'plotMean',false,...
        'faceColor',[.5 .75 1],'alpha',1,'lineColor',[.3 .5 1]); 
hold on;
med = mean( double(vdists),1 );
plot(1:4,med,'k.','MarkerSize',30);
ylim([0,70]);
ylabel('transcription rate');
xlabel('p1 WT   p2 WT   p1 Del   p2 Del');
set(gcf,'color','w');

%%
% SaveFigure(f4,'name','fig6_P1vsP2_violin','formats',{'png','eps'});