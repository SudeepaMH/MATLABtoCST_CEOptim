function [X_opt, y_opt, time, flag] = soce(CST, l, u, options, g, gamma)
%SOCE Cross-entropy single-objective optimization
%   Applies the cross-entropy method to a single-objective optimization 
%   problem
%
%   Syntax:
%      [X_opt, y_opt, time, flag] = SOCE(f, l, u, options, g, gamma)
%
%   Arguments:
%      f [function reference or inline function] Objective functions. It
%         must return a scalar value.
%      l [real, n-values vector] Lower bound of the decision variables.
%      u [real, n-values vector] Upper bound of the decision variables.
%      options [structure] Optimization options. It can be defined by using
%         the function CEOPTDEF.
%      g [function reference or inline function]. Constraints functions. It
%         must return a p-values vector. For unconstrained problems this
%         argument should be empty.
%      gamma [real, p-values vector] Constraints penalty factors. For
%         unconstrained problems this argument should be empty.
%
%   Outputs:
%      x_opt [real, n-values vector] Values of the decision variables at the
%         optimal point
%      y_opt [real, scalar] Value of the objective function at the optimal point.
%      time [real, scalar] Execution time
%      flag [string] Text describing why the program stopped.
% 
%   Ref.: Rubinstein, R.Y.; Kroese, D.P., 2004. "A Tutorial Introduction 
%         to the Cross-Entropy Method”. The Cross-Entropy Method. New York:
%         Springer, ISBN 978-1-4419-1940-3, pp. 29-58.
%
%   Cross-Entropy Optimisation Toolbox
%   (c) 2018, GAMHE <http://www.gamhe.eu> and CEFAS <http://cefas.umcc.cu>
%

   % Verifying the arguments values (CST)
    if (nargin < 1), error('Not enough input arguments.'); end
    if isempty(CST), error('Argument ''CST'' should not be empty.'); end
    %    if ~(isa(f, 'function_handle') || isa(f, 'inline')), error('Argument ''f'' should be a function reference or an inline funtion.'); end
   % L
   if (nargin < 2), error('Not enough input arguments.'); end
   if isempty(l), error('Argument ''l'' should not be empty.'); end
   if ~isnumeric(l), error('Argument ''l'' should be a numeric value.'); end
   if (size(l, 1) > 1), error('Argument ''l'' should be a single row vector.'); end
   % U
   if (nargin < 3), error('Not enough input arguments.'); end
   if isempty(u), error('Argument ''u'' should not be empty.'); end
   if ~isnumeric(u), error('Argument ''u'' should be a numeric value.'); end
   if (size(u, 1) > 1), error('Argument ''u'' should be a single row vector.'); end
   if ~(all(size(l) == size(u))), error('Arguments ''l'' and ''u'' should be vectors with the same length.'); end
   if any(u <= l), error('All values of argument ''u'' must be higher than the corresponding values of ''l''.'); end
   % OPTIONS
   if (nargin < 4), options = ceoptdef(); end
   if isempty(options), error('Argument ''options'' should not be empty.'); end
   if ~isstruct(options), error('Argument ''options'' should be a structure.'); end
   includedFields = { 'Generations', 'PopulationSize', 'EliteRatio', ...
                      'Display', 'ConvergenceLimit', 'Alpha', 'Beta', ...
                      'ExpFactor', 'SmoothingType' };
   for i = 1 : length(includedFields)
      if ~isfield(options, includedFields{i}), error(['Argument ''options'' should contain a field called ''', includedFields{i}, '''.']); end
   end
   % G
   if (nargin < 5), g = []; end
   if ~(isempty(g) || isa(g, 'function_handle') || isa(g, 'inline')), error('Argument ''g'' should be a function reference or an inline funtion.'); end
   % GAMMA
   if (nargin < 6), gamma = []; end
   if (isempty(gamma) && ~isempty(g)), error('Argument ''gamma'' should not be empty if argument G is not.'); end
   if ~isempty(g)
      if ~isnumeric(gamma), error('Argument ''gamma'' should be a numeric value.'); end
      if (size(g, 1) > 1), error('Argument ''gamma'' should be a single row vector.'); end
   end

   % Computing the number of decision variables and constraints
   n = length(l);
   p = length(gamma);
   
   % Loading options
   N = options.Generations;
   Z = options.PopulationSize;
   rho = options.EliteRatio;
   Ne = ceil(rho.*Z);
   display = options.Display;
   epsilon_max = options.ConvergenceLimit;
   alpha = options.Alpha;
   beta = options.Beta;
   q = options.ExpFactor;
   smoothingType = options.SmoothingType;

   % Initialization
   sigma(1 : n) = 2.*(u - l);
   mu(1 : n) = l + rand(1, n).*(u - l);
   mu_last = mu;
   sigma_last = sigma;
   X_opt = zeros(1, n);
   t = 0;
   Q = zeros(Z, n + 1);
   L = [];
   time = now;
   %% Debug handles
   isdebug = 0;
   %%
   optimizing_freq_index = [1,1,1,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8];
   optimizing_freq = [1563,1575,1587,1472,2310,2335,2360,2400,2441.75,2483.5,2400,2442,2484,5150,5250,5350,5725,5787.5,5850,5725,5800,5875];
   
   %%
   target_score = Z+length(optimizing_freq);
    %%
    if (strcmpi(display, 'iter'))
       fprintf('Generation  Evaluations       Time [s]           Convergence     Best f(x)     Mean f(x)\n');
    end
   
    %% Open text file for logging
    fid = fopen('generation_logger.txt','wt');
    %fclose(fid);
    
   % Iterations
   while 1
      % Increasing the epoch counter
      t = t + 1;
      % Updating means and standard deviations
      mu = alpha.*mu + (1 - alpha).*mu_last;
      beta_mod = beta - beta.*(1 - 1./t).^q;
      if (smoothingType == 1)
         sigma = beta_mod*sigma + (1 - beta_mod).*sigma_last;
      else
         sigma = alpha.*sigma + (1 - alpha).*sigma_last;
      end
      % Computing and recording variables values in the new population
      Q(:, 1 : n) = randnt(Z, n, mu, sigma, l, u); 
      
      fprintf(fid, '===== Generation %d =====\n', t);
      fprintf(fid, 'mu = %.8f\n', mu);
      fprintf(fid, 'sigma = %.8f\n', sigma);
      fprintf(fid, '**** Simulation Results ****\n');
      fprintf(fid, 'diople19_spacing \t diople18_spacing \t diople17_spacing \t diople16_spacing \t diople15_spacing \t diople14_spacing \t diople13_spacing \t d12_spacing \t d11_spacing \t Resonance Sum \t Total Length\n');        
      for i=1:Z
            d19_spacing_this = Q(i,1);
            d18_spacing_this = Q(i,2);
            d17_spacing_this = Q(i,3);
            d16_spacing_this = Q(i,4);
            d15_spacing_this = Q(i,5);
            d14_spacing_this = Q(i,6);
            d13_spacing_this = Q(i,7);
            d12_spacing_this = Q(i,8);
            d11_spacing_this = Q(i,9);
            
            if(isdebug == 0)
                CST.changeParameterValue('dipole19_spacing',d19_spacing_this);
                CST.changeParameterValue('dipole18_spacing',d18_spacing_this);
                CST.changeParameterValue('dipole17_spacing',d17_spacing_this);
                CST.changeParameterValue('dipole16_spacing',d16_spacing_this);
                CST.changeParameterValue('dipole15_spacing',d15_spacing_this);
                CST.changeParameterValue('dipole14_spacing',d14_spacing_this);
                CST.changeParameterValue('dipole13_spacing',d13_spacing_this);
                CST.changeParameterValue('dipole12_spacing',d12_spacing_this);
                CST.changeParameterValue('dipole11_spacing',d11_spacing_this);
            end
            if(isdebug == 0)   
                CST.runSimulation;           
                [freq,s11,type] = CST.getSParameters('S11');
                S11_mag = abs(s11);
                s11_dB = 20*log10(S11_mag);
                s11_plot = [freq; s11_dB];
                s11AtOptimizingFreq = [];
                resonanceSum = 0;
                for j = 1 : length(optimizing_freq)
                    isResonate = 0;
                    thisOptimizingFreq = optimizing_freq(j) ;
                    [~,indexofClosestVal] = (min(abs(freq - thisOptimizingFreq)));
                    if(s11_dB(indexofClosestVal) < -10.0)
                        isResonate = 1;
                    end
                    resonanceSum = resonanceSum + isResonate;
                    thisSearch = [optimizing_freq_index(j) thisOptimizingFreq s11_dB(indexofClosestVal) isResonate]';
                    s11AtOptimizingFreq = [s11AtOptimizingFreq thisSearch];
                end
                Q(i, n + 1) = resonanceSum;
                thisAntennaLength = getParameterValue(CST,'substrate_length');
                L = [L; i thisAntennaLength];
            end
            if(isdebug == 1)
                Q(:, n + 1) =  1 + (10-1).*rand(1);
                %Q(:, n + 1) =   1 + (10-1).*rand(n,1);
                %Q(:, n + 1) = f(Q(:, 1 : n)); 
            end
            fprintf(fid, '%.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f\n', d19_spacing_this, d18_spacing_this, d17_spacing_this, d16_spacing_this, d15_spacing_this, d14_spacing_this, d13_spacing_this, d12_spacing_this, d11_spacing_this, Q(i, n + 1), L(i,2));
            i=i+1;                  
      end
      %%  Sort the antenna by length and compute the final score for each antenna.
      % Sorting the antenna by length. 
      % Shortest Antenna gets 20 points, Longest Antenna gets 1 point
      sortedL = sortrows(L, 2,'ascend');
      lengthScores = [Z:-1:1]';
      sortedL = [sortedL lengthScores];
      sortedL = sortrows(sortedL, 1,'ascend');
      
      % Computing final score of the antenna [LengthScore + Resonance Sum]
      Q(:,n+1) = Q(:,n+1)+sortedL(:,3)
      %% Computing and recording constraints values in the new population
      if (p > 0)
         G = g(Q(:, 1 : n));         
         Q(:, n + 1) = Q(:, n + 1) + sum(max(G, 0)*diag(gamma), 2);
      end
      % Sorting the solutions by their performance
      Q = sortrows(Q, n + 1,'descend');
      % Updating the best solution
      if (Q(1, n + 1) > target_score)
         y_opt = Q(1, n + 1);
         X_opt = Q(1, 1 : n);
      end
      % Extractiing the elitist population
      E = Q(1 : Ne, :);
      fprintf(fid, '**** Elite Population ****\n');
      fprintf(fid, 'diople19_spacing \t diople18_spacing \t diople17_spacing \t diople16_spacing \t diople15_spacing \t diople14_spacing \t diople13_spacing \t d12_spacing \t d11_spacing \t Final Score\n');        
      fprintf(fid,'%.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \t %.8f \n',E.');
      % Considering the stop conditions
      mu_last = mu;
      sigma_last = sigma;
      % Showing process progress
      if (strcmpi(display, 'iter'))
         fprintf(' %9d   %10d  %13.4f     %8.2e/%8.2e  %12.4e  %12.4e\n', ...
                 t, t*Z, (now - time)*24*60*60, max(sigma), epsilon_max, ...
                 max(E(:, n + 1)), mean(E(:, n + 1)));
      end      
      if (t >= N)
         % Stop if raching the maximum number of epochs
         time = (now - time)*24*60*60;
         flag = 'Maximum epoch reached';
         break;
      elseif (max(sigma) <= epsilon_max)
         % Stop if all the standar deviations are lower than epsilon_max
         time = (now - time)*24*60*60;
         flag = 'Convergence limit reached';
         break;
      end
      % Computing the new means and standard deviations
      mu = mean(E(:, 1 : n));
      sigma = std(E(:, 1 : n));
      fprintf(fid, 'mu elite = %.8f\n', mu);
      fprintf(fid, 'sigma elite = %.8f\n', sigma);
      fprintf(fid, '\n');
      sortedL = [];
      L = [];
   end
  
   % Showing process completion
   if (strcmpi(display, 'iter') ||strcmpi(display, 'final'))
     fprintf('Optimization process completed.\n');
     fprintf('   Stop due to: %s.\n', flag);
     fprintf('   Epoch number: %d.\n', t);
%      fprintf('   Objective function f(x) = %f.\n', y_opt);
%      for i = 1 : n
%          fprintf('   Decision variable x%d = %f.\n', i, X_opt(i));
%      end
%      if (p > 0)
%          G_opt = g(X_opt);
%          for i = 1 : p
%              fprintf('   Constraint g%d(x) = %f.\n', i, G_opt(i));
%          end
%      end
   end 
end

