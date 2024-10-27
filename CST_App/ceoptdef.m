 function options = ceoptdef(varargin)
%CEOPTDEF Sets the options for the cross-entropy method
%   Creates a structure containing the default values of the options for
%   the cross-entropy method.
%
%   Syntax:
%      options = CEOPTDEF()
%         Creates an OPTIONS structure with the default settings.
%      options = CEOPTDEF('param1', value1, 'param2', value2...)
%         Creates an OPTIONS structure with the PARAM-VALUES pairs
%         setting.
%      options = CEOPTDEF(options, 'param1', value1, 'param2', value2...)
%         Modify the argument OPTIONS by changing the PARAM-VALUES pairs
%         setting.
%
%   Cross-entropy parameters {default values}
%   Single and multi-objective optimization
%      Generations [integer scalar {100}] The maximum number of epoch
%            after which the optimization process stops.
%      PopulationSize [integer scalar {50}] The number of solutions
%            simultaneously evaluated.    
%      EliteRatio [real scalar {0.1}] Ration between the elitist population
%            and the working of population.
%      Display ['off' | 'iter' | 'diagnose' | {'final'}] Level of display.
%   Only single-objective optimization
%      ConvergenceLimit [real scalar {0.1}] Maximum change in the
%            standard deviations.
%      Alpha [real scalar (0...1) {0.7}] Smoothing factor of means and
%            standard deviations.
%      Beta [real scalar (0...1) {0.7}] Factor of dynamical smoothing of
%            means and standard deviations.
%      ExpFactor [real scalar {50}] Exponent of dynamical smoothing of
%            means and standard deviations.
%      SmoothingType [{'dynamical'} | 'fixed'] Smotthing type.
%   Only multi-objective optimization
%      HistogramIntervals [integer scalar {5}] Number of intervals in the
%            objective-space frequency histogram.
%
%   Cross-Entropy Optimisation Toolbox
%   (c) 2018, GAMHE <http://www.gamhe.eu> and CEFAS <http://cefas.umcc.cu>
%

   % Default fields names and values
   defFieldNames{ 1} =  'Generations';         defFieldValues{ 1} = 100;
   defFieldNames{ 2} =  'PopulationSize';      defFieldValues{ 2} = 50;
   defFieldNames{ 3} =  'EliteRatio';          defFieldValues{ 3} = 0.1;
   defFieldNames{ 4} =  'Display';             defFieldValues{ 4} = 'final';   
   defFieldNames{ 5} =  'ConvergenceLimit';    defFieldValues{ 5} = 0.1;
   defFieldNames{ 6} =  'Alpha';               defFieldValues{ 6} = 0.8;
   defFieldNames{ 7} =  'Beta';                defFieldValues{ 7} = 0.7;
   defFieldNames{ 8} =  'ExpFactor';           defFieldValues{ 8} = 50;
   defFieldNames{ 9} =  'SmoothingType';       defFieldValues{ 9} = 'dynamical';
   defFieldNames{10} =  'HistogramIntervals';  defFieldValues{10} = 5;

   
   % Checking if there is not input arguments
   if (nargin == 0)
      options = struct(); 
      for i = 1 : length(defFieldNames)
         options = setfield(options, defFieldNames{i}, defFieldValues{i});
      end
      return;
   end
   
   % Checking if the first argument is a structure
   if isstruct(varargin{1})
      n = 1;
      options = varargin{1};
      for i = 1 : length(defFieldNames)
         if ~isfield(options, defFieldNames{i})
            options = setfield(options, defFieldNames{i}, defFieldValues{i});
         end
      end
   else
      n = 0;
      options = struct(); 
      for i = 1 : length(defFieldNames)
         options = setfield(options, defFieldNames{i}, defFieldValues{i});
      end
   end
   
   % Checking the name-value pairs structure of the arguments
   if rem(nargin + n, 2) ~= 0
      error('Arguments must occur in name-value pairs.');
   end
   
   % Getting options
   for i = n + 1 : nargin
      if rem(i + n, 2)
         if ~ischar(varargin{i})
            error('Expected argument %d to be a string parameter name.', i);
         end
         if ~isfield(options, varargin{i})
            error('''%s'' is not a valid option name.', varargin{i});
         end
         options = setfield(options, varargin{i}, varargin{i + 1});
      end
   end  
end
