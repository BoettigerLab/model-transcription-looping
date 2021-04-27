% Stochastic PTM version of the model.  Produces the figures for the PTM
% violin plots in Figure 3. 
% 
%
% Corresponding ODE model:
% da/dt = (Et*(k1*k2*r2*a.^2 + k1*r1*tau2*a + r0*tau1*tau2))./(k1*k2*a.^2 + k1*tau2*a + tau1*tau2) +loopRate;

% Model Parameters
k1 = .003; % Rate of conversion from state 0 to state 1.
k2 = .5; % Rate of conversion from state 1 to state 2. 
tau1 = 20; % Rate of conversion from state 1 to state 0.
tau2 = .05; % Rate of conversion from state 2 to state 1.
hdac = 2; % Rate of tag removal. 
r0 = .01; % Basal tagging rate (state 0)
r1 = .5; % Stimulated tagging rate (state 1) 
r2 =5.0; % Maximally stimulated tagging rate (state 2)
loopRate = 9; % .65 ;% .95; % the linear shift to the sigmoid. % 20% change, 3.5 to 4.5  
Et= 20; % Maximum tag-transferases associated with promoter

% initial conditions
ac_init = 0; %  No tags (40 = start high)
Pr_init = 0; % all enzymes in the unstimulated state 0

% E-P loop rates explored
loopRates = 0:.2:20; 

% number of simulations, length of simulation, time scale of sampling
tSteps = 3e5; % length of simulation
nRuns = 3e3; % number of replicates
t_scale = 1000; % time scale of sampling
numLoops = length(loopRates);

%% Run the simulation
saveFolder = SetFigureSavePath('U:\Alistair\Analysis\2020-04-23_JordanLoopSims\');
saveName = [saveFolder,'loopRatesVsExpStoch.mat'];
% only need to repeat the simulation if a saved file is not found
if exist(saveName,'file') == false    
    acVals2 = cell(numLoops,1);
    for r=1:numLoops
        disp(r/numLoops);
        loopRate = loopRates(r);  
        acR = [r0 r1 r2];  % acetylation rates for the different enhancer states 
        compTx = false;

        nZymes =Et;
        clear ac_values;
        ac_values = zeros(nRuns,tSteps,'uint8');

        parfor n=1:nRuns
            acR = [r0 r1 r2];
            ac = ac_init; % 0 ... inf
            Pr = Pr_init*ones(1,nZymes);  % 1, 2, 3
            for t=1:tSteps % straight discreet time simulation

                % enhancer mediated actyl
                if rand*t_scale < loopRate
                   ac = ac + 1; 
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
            ac_values(n,t) = ac;
            end
        end
        acVals2{r} = ac_values;
        acv = ac_values(:,end-500:end);
    end

%% Save or load results
    save(saveName,'acVals2','expR','-v7.3');
else
    load([saveFolder,'loopRatesVsExpStoch.mat'],'acVals2');
end
%% plot stoch results

% reorganizing data into matrix for plotting
numLoops = length(acVals2);
[nRuns, tSteps] = size(acVals2{1});
lowDist = zeros(numLoops,nRuns); 
for r=1:numLoops
    lowDist(r,:) = acVals2{r}(:,end);
end

% plot results as a violin plot
f3 = figure(3); clf;
subsample = 1:5:numLoops-10;
violin(lowDist(subsample,:)','bandwidth',.5,'plotMean',false,'faceColor',[.5 .5 1],'alpha',1); hold on;
lowMed = median(lowDist(subsample,:),2);
plot(lowMed,'b'); hold on;
ylim([0,70]);
set(gca,'xTickLabel',loopRates(subsample));
set(gcf,'color','w');
xlabel('loop rate');
ylabel('transcription rate');
%%
% SaveFigure(f3,'name','violin_stochBistable_low','formats',{'png','eps'});


