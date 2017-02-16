classdef GA2 < handle
   % Genetic algorithm
   %
   % USAGE:
   %     obj = GA2(popSize, ngenes, objFun, OPTIONAL objFunParam, OPTIONAL bits)
   %     obj.optimize( maxGens );
   %     obj.step();
   % ARGS
   %     popSize     - population size 
   %     ngenes      - number of parameters of the objective function
   %     objFun      - handle to the objective function of the form
   %                   `fitness = objFun( parameters, objFunParam)` 
   %     objFunParam - OPTIONAL - parameters to pass to objective function
   %     bits        - bit resolution of the parameters - OPTIONAL, defaults to 64 bits
   % RETURNS
   %     GA2 object with these PROPS (selection):
   %        eliteIndivParam - current best set of parameters
   %        avgFitnessHist  - avg. fitness of the solutions in the pop for
   %                          each generation
   %        maxFitnessHist  - avg. fitness of the solutions in the pop for
   %                          each generation
   
   % modified from: 
   %     [TurboGA](http://de.mathworks.com/matlabcentral/fileexchange/24053-turboga--a-simple-genetic-algorithm-with-a-powerful-performance-enhancing-tweak) 
   % last mod 20151228 janc
   
   % IMPLEMENTATION DETAILS
   % _genome_ - 64bit binary vectors
   % _clamping_ -
   % _crossover_ - no crossover/1pt crossover/uniform crossover (default)
   % _point mutations_ - pretty standard
   % _selection_: Stochastic Universal Sampling (default) or Roulette Wheel Sampling
   % _sigma scaling_: selection based on normalized (zero-mean, unit-variance) fitness values
   % _elite pool_: keep the n best solutions in a separate pool
   
   properties
      % clamping
      clampingFlag = 1;% Use a mechanism called clamping (see http://www.cs.brandeis.edu/~kekib/GAWorkings.pdf for details)
      flagFreq = 0.01;
      unflagFreq = 0.1;
      flagPeriod = 200;
      flaggedGens
      
      % crossover
      crossoverType = 2;      % 0 - no crossover, 1 - 1pt crossover, 2 - uniform crossover, if clampingFlag = 1, crossoverType should be 2
      probCrossover = .96;    % The probability of crossing over.
      
      % point mutations
      probMutation = 0.003;   % The mutation probability (per bit). If clampingFlag = 1, probMutation should not be dependent on len, the length of the genomes.
      
      % selection
      SUSFlag = 1;            % Use Stochastic Universal Sampling (pg 168 of M. Mitchell's GA book), almost always improves performance
      sigmaScalingFlag = 1;   % Sigma Scaling is described on pg 168 of M. Mitchell's GA book. It often improves GA performance.
      sigmaScalingCoeff = 1;  % Higher values = > less fitness pressure
      sigma                   %
      
      % elite pool
      elitePoolFlag = 1;   % never loose the best solution by always keeping (and never mutating) the ELITEPOOLSIZE best solutions
      elitePoolSize        % number of best individuals to clone to the next generation
      elitePool            % the gnomes of the `elitePoolSize` best inviduals
      eliteFitness = -realmax; %init to smallest possible value
      eliteIndivGenome;    % binary genome of the best individual
      eliteIndivParam;     % parameters of the best individual
      
      avgFitnessHist, maxFitnessHist % vectors for recording the average and maximum fitness in each generation
      
      verboseFlag = 1;% display details of each generation; 0: run quietly
      
      bits = 64;  % bit resolution of the genome
      fitnessVals % holds the current generations fitness values
      popSize     % population size
      ngenes      % number of parameters
      len         % genome length = bits*ngenes
      objFun      % handle to the objective function - CONDITIONS
      objFunParam % parameters
      pop         % genomes of population
      gen         % generation counted
   end
   
   
   methods
      function obj = GA2(popSize, ngenes, objFun, objFunParam, varargin)
         
         obj.popSize = popSize;
         obj.ngenes = ngenes;
         
         % optional: bit resolution of the parameters
         if ~isempty(varargin)
            obj.bits = varargin{1};
         end
         % optional: parameters to pass to objective function
         if ~exist('objFunParam', 'var')
            objFunParam = [];
         end
         
         % objective function handle and parameters
         obj.objFun = objFun;
         obj.objFunParam = objFunParam;
         
         % genome size
         obj.len = obj.bits*obj.ngenes;
         obj.flaggedGens = -ones(1,obj.len);
         
         obj.gen = 0;
         % the population is a popSize by len matrix of randomly generated binary values
         obj.pop = rand(obj.popSize,obj.len)<.5;
         
         % number of best individuals to clone to the next generation
         obj.elitePoolSize = limit(ceil(round(obj.popSize/50)/2), 1, inf);
         
         % preallocate vectors for recording the average and maximum fitness in each generation
         obj.avgFitnessHist = zeros(0);
         obj.maxFitnessHist = zeros(0);
      end
      
      function optimize(obj, nSteps)
         % obj.optimize([nSteps=100])
         if nargin == 1
            nSteps = 100;
         end
         
         for stp = 1:nSteps
            obj.step();
         end
         
      end
      
      function step(obj)
         % obj.step()
         % perform single optimization step
         if obj.verboseFlag
            tic
         end
         % get the parameters from the binary genome
         param = binVec2decVec(obj.pop, obj.bits);
         % evaluate the objective function
         obj.fitnessVals = feval(obj.objFun, param, obj.objFunParam);
         
         % fix invalid solutions
         obj.fitnessVals(isinf(obj.fitnessVals)) = nan;%realmin('double');
         obj.fitnessVals(isnan(obj.fitnessVals)) = nanmin(obj.fitnessVals);
         
         obj.avgFitnessHist(1,obj.gen+1) = nansum(obj.fitnessVals./obj.popSize);
         [val, idx] = sort(obj.fitnessVals);
         obj.maxFitnessHist(1,obj.gen+1) = val(end);
         % keep the N best individuals in a separate pool
         maxIndex = idx(end+1-obj.elitePoolSize:end);
         obj.elitePool = obj.pop(maxIndex,:);
         % if we found a new best individual, add to elite
         if obj.eliteFitness < obj.maxFitnessHist(1,obj.gen+1)
            obj.eliteFitness = obj.maxFitnessHist(1,obj.gen+1);
            obj.eliteIndivGenome = obj.pop(idx(end),:);
            obj.eliteIndivParam = param(idx(end),:);
         end
            
         % sigma scaling
         if obj.sigmaScalingFlag
            obj.sigma = obj.popSize*nanstd(obj.fitnessVals./obj.popSize);
            if obj.sigma~=0
               obj.fitnessVals = 1 + (obj.fitnessVals-nansum(obj.fitnessVals/obj.popSize))/...
                  (obj.sigmaScalingCoeff*obj.sigma);
               obj.fitnessVals(obj.fitnessVals <= 0) = 0;
            else % fallback
               obj.fitnessVals = ones(1, obj.popSize);
            end
         end
         
         % Normalize the fitness values to unit-sum and then create an array with the
         % cumulative normalized fitness values (the last value in this array will be 1)
         cumNormFitnessVals = cumsum(obj.fitnessVals(~isnan(obj.fitnessVals))/(nansum(obj.fitnessVals)+eps));
         
         % determine which genes to clamp
         bitFreqs = sum(obj.pop)/obj.popSize;
         if obj.clampingFlag
            lociToFlag = abs(0.5-bitFreqs)>(0.5-obj.flagFreq) & obj.flaggedGens<0;
            obj.flaggedGens(lociToFlag) = 0;
            lociToUnflag = abs(0.5-bitFreqs)<0.5-obj.unflagFreq ;
            obj.flaggedGens(lociToUnflag) = -1;
            flaggedLoci = obj.flaggedGens>= 0;
            obj.flaggedGens(flaggedLoci) = obj.flaggedGens(flaggedLoci)+1;
            mutateLocus = obj.flaggedGens<= obj.flagPeriod;
         end
         
         % Use fitness proportional selection with Stochastic Universal or Roulette
         % Wheel Sampling to determine the indices of the parents of all crossover operations
         if obj.SUSFlag
            markers = rand(1,1)+(1:obj.popSize)/obj.popSize;
            markers(markers>1) = markers(markers>1)-1;
         else
            markers = rand(1,obj.popSize);
         end
         [~, parentIndices] = histc(markers,[0 cumNormFitnessVals]);
         parentIndices = parentIndices(randperm(obj.popSize));
         
         % determine the first parents of each mating pair
         firstParents = obj.pop(parentIndices(1:obj.popSize/2),:);
         % determine the second parents of each mating pair
         secondParents = obj.pop(parentIndices(obj.popSize/2+1:end),:);
         
         % create crossover masks
         if obj.crossoverType == 0
            masks = false(obj.popSize, obj.len);
         elseif obj.crossoverType == 1
            masks = false(obj.popSize/2, obj.len);
            temp = ceil(rand((obj.popSize)/2,1)*(obj.len-1));
            for i = 1:(obj.popSize)/2
               masks(i,1:temp(i)) = true;
            end
         else
            masks = rand((obj.popSize)/2, obj.len)<.5;
         end
         
         % determine which parent pairs to leave uncrossed
         reprodIndices = rand((obj.popSize)/2,1)<1-obj.probCrossover;
         masks(reprodIndices,:) = false;
         
         % implement crossover
         firstKids = firstParents;
         firstKids(masks) = secondParents(masks);
         secondKids = secondParents;
         secondKids(masks) = firstParents(masks);
         obj.pop = [firstKids; secondKids];
         
         % point mutations
         masks = rand(obj.popSize, obj.len)<obj.probMutation;
         if obj.clampingFlag
            masks(:,~mutateLocus) = false;
         end
         obj.pop = xor(obj.pop,masks);
         
         % always keep the elite in there...
         obj.pop(maxIndex,:) = obj.elitePool; 
         
         % update internal generation counter
         obj.gen = obj.gen+1;
         
         % display some info
         if obj.verboseFlag
            fprintf('      gen = %.3d:  avgFitness = %3.3f  maxFitness = %3.3f  took %3.2fs.\n',...
               obj.gen, obj.avgFitnessHist(1, obj.gen), obj.maxFitnessHist(1, obj.gen), toc)
         end
      end
      
   end % methods
end % classdef
