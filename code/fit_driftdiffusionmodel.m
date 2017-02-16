cc()
fprintf('   loading data.\n')
d = load('data.mat');
%%
pa.xvrun = 1; % stimulus index to hold out for cross-validation - iterate this over all 218 stimuli 
pa.envSyllType = d.stim;  % subunit sequences for each stimulus
pa.meanResp = d.meanResp; % behavioral responses
pa.stimType = d.stimType; % subunit pair tested - for plotting purposes 
pa.stis = length(pa.meanResp);
clear d;

% upper and lower bounds for the parameters - helps speed up fitting - the
% results are qualitatively robust to the exact choise of parameter bounds
pa.paramLabel = {'w-dorsatus','w-gap','w-mollis','w-onset8dB','w-pause4ms','w-onset3dB','w-onset9dB', 'w-accentOffset','w-accentOnset', '\sigma', '\theta_+', '\theta_-'};
pa.lb = [-100 -100 -100 -100 -100 -100 0 -100 0 0 0 0];
pa.ub = [0 0 0 0 0 0 50 0 50 200 500 1000];

% function that generates the training error for each individual
% this implements the actual drift-diffusion model
pa.objFunInd = @LEI2_perfectMulti; 
pa.batch = 1.0; % fraction of data to use for training in each iterations

% generate integrator noise - this speeds up training
pa.noiseRuns = 2000; % noise repeats
pa.cumNoiseSize = [33, pa.noiseRuns,pa.stis];
pa.cumNoise = cumsum(randn( pa.cumNoiseSize), 1);

paTrain = pa; % copy for 
fprintf('   holding out stimulus %d for later crossvalidation.\n', pa.xvrun);
% hold-out test stimulus
paTrain.stis = paTrain.stis-1;
paTrain.meanResp(pa.xvrun) = [];
paTrain.envSyllType(:,pa.xvrun) = [];

fprintf('   initializing genetic algorithm.\n');
popSize = 200; % number of individuals in the population of the genetic algorithm
maxGens = 500; % number of generations to fit
objFun = @LEIpop;% the function that iterates over all individuals
ngenes = length(pa.paramLabel); % number of genes (=parameters)

pop = GA2(popSize, ngenes, objFun, paTrain); % init the GA fitter - see comments in GA2.m
%%
fprintf('   fitting model for %d generations - this will take a long time...\n', maxGens)
for gen = 1:maxGens
   pop.step(); % run one generation - this is doing the actual work
   
   % process and plot results:
   % rescale best solution form normalized units (-1 1) to (lb ub)
   for i = 1:size(pop.eliteIndivParam,2)
      best.param(:,i) = (pop.eliteIndivParam(:,i)+1)/2*(pa.ub(i) - pa.lb(i)) + pa.lb(i);
   end
   % get prediction for current best solution
   [best.er, best.pred, best.intEvidence] = pa.objFunInd(best.param, pa);

   % plot intermediate results
   subplot(231)
   plot(pop.maxFitnessHist)
   hline(max(pop.maxFitnessHist))
   ylabel('mse')
   xlabel('generation')
   title('history')
   set(gca, 'XTick', 1:gen)

   subplot(234)
   bar(pop.eliteIndivParam)
   axis('tight')
   hline(0)
   xlabel('parameter #')
   ylabel('values in normalized units')
   set(gca, 'YLim', [-1 1], 'XTick', 1:length(pa.paramLabel), ...
            'XTickLabel', pa.paramLabel, 'XTickLabelRotation', 90)

   subplot(1,3,2:3)
   gscatter(best.pred, pa.meanResp, pa.stimType, [],'o')
   title(sprintf('generation %d, r^2=%1.2f', gen, rsq(best.pred, pa.meanResp)))
   set(gca, 'XLim', [0 1], 'YLim', [0 1])
   axis('square')
   dline()
   xlabel('prediction')
   xlabel('behavior')
   drawnow
end

% display final performance metrics
disp(rsq(best.pred, pa.meanResp))
disp(1-best.er/std(pa.meanResp))

% remove cumNoise to save space/bandwidth
pa.cumNoise = [];
pop.objFunParam.cumNoise = [];

% save results
% pa - training data and parameters
% pop - the GA object
% best - current best solution
save(sprintf('xv_%03d.mat', pa.xvrun), 'pa', 'pop', 'best')
