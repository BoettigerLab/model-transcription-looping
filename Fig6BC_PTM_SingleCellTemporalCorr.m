
%% Fig 5B
saveFolder = SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Images\');
% save([saveFolder,'oddsRatioEmbryoData.mat'],'values','errors','contactFreq');

load([saveFolder,'oddsRatioEmbryoData.mat'],'values','errors','contactFreq');
[c,i] = sort(contactFreq);

% [c,i] = sort(values);
f1 = figure(1); clf;
ploterr(c,values(i),[],{values(i)-errors(i),values(i)+errors(i)},'.','color','k');
hold on; plot([0,100],[1,1],'k--');
ylim([0,2]);
xlim([.12,.38]);
xlabel('E-P contact frequency');
ylabel('Odds Ratio');
title('Odds of transcription given contact');


%% Fig 5C
load([saveFolder,'loopResultsStartLow.mat'],'loopVals','acVals','loopRates');
% plot
orON = zeros(20,1);
orciON = zeros(20,2);
or_stdev = zeros(20,1);
for r=1:20 % r = 7;
    lowCell_loops = double(cat(2,loopVals{r/.2}{:}));
    lowCell_expr = double(cat(2,acVals{r/.2}{:}));
    lowCell_loops = lowCell_loops(1000:2000,1:100);
    lowCell_expr = lowCell_expr(1000:2000,1:100);
    lowON = false(size(lowCell_loops));
    lowFire = lowON;
    lowON(lowCell_expr > 20) = true; % 20 lowON(lowCell_expr > 50) = true;
    lowFire(lowCell_loops > 1) = true; %0, 1   4 lowFire(lowCell_loops > 10) = true;
    [orON(r),orciON(r,:),or_stdev(r)] = OddsRatioCI(lowFire(:),lowON(:),'iters',10+ceil(60/r));
end

f2 =figure(2); clf; set(gcf,'color','w');
ploterr((1:20)',orON,[],{orON-or_stdev,orON+or_stdev},'.','color','k');
xlabel('E-P contact frequency');
ylabel('Odds Ratio');
title('Odds of transcription given contact');
hold on; plot([0,100],[1,1],'k--');
xlim([3.5,16.5]); ylim([0,2]);

%% Fig 5D
