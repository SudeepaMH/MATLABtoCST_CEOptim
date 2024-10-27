function Z = randnt(n, m, mu, sigma, l, u)
%RANDNT Normally distributed truncated pseudorandom numbers
%   Generates an N-by-M matrix of numbers drawn form the normal
%   distributions with means MU and standard deviations SIGMA and truncated 
%   to the intervals defined by the lower bounds L and the upper bounds U.
%
%   Syntax:
%      Z = RANDNT(n, m, mu, sigma, l, u)
%
%   Arguments:
%      n [integer, scalar] Number of rows in the generated matrix.
%      m [integer, scalar] Number of columns in the generated matrix.
%           Each column corresponds to one normal truncated distribution.
%      mu [real, M-elements row vector] Mean value of each normal
%           truncated distribution (default MU = 0).
%      sigma [real, M-elements row vector] Standard deviation of each
%           normal truncated distribution (default SIGMA = 1).
%      l [real, M-elements row vector] Lower bound of each normal
%           truncated distribution (default L = -Inf).
%      u [real, M-elements row vector] Upper bound of each normal
%           truncated distribution (default U = Inf).
%
%   Outputs:
%      Z [real, N-by-M matrix] Normally distributed truncated numbers.
%
%   Cross-Entropy Optimisation Toolbox
%   (c) 2018, GAMHE <http://www.gamhe.eu> and CEFAS <http://cefas.umcc.cu>
%

   % Setting default values for the arguments
   if (nargin < 1), n = 1; end
   if (nargin < 2), m = 1; end
   if (nargin < 3), mu = zeros(1, m); end
   if (nargin < 4), sigma = ones(1, m); end
   if (nargin < 5), l = -Inf.*ones(1, m); end
   if (nargin < 6), u = Inf.*ones(1, m); end
   
   % Verifying the arguments values (N)
   if isempty(n), error('Argument ''n'' should not be empty.'); end
   if ~isnumeric(n), error('Argument ''n'' should be a numeric value.'); end
   if (mod(n, 1) > 0), error('Argument ''n'' should be an integer value.'); end
   if (max(size(n)) > 1), error('Argument ''n'' should be a scalar value.'); end
   if (n < 1), error('Argument ''n'' should be a pisitive value.'); end
   % M
   if isempty(m), error('Argument ''m'' should not be empty.'); end
   if ~isnumeric(m), error('Argument ''m'' should be a numeric value.'); end
   if (mod(m, 1) > 0), error('Argument ''m'' should be an integer value.'); end
   if (max(size(m)) > 1), error('Argument ''m'' should be a scalar value.'); end
   if (m < 1), error('Argument ''m'' should be a positive value.'); end
   % MU
   if isempty(mu), error('Argument ''mu'' should not be empty.'); end
   if ~isnumeric(mu), error('Argument ''mu'' should be a numeric value.'); end
   if (size(mu, 1) > 1), error('Argument ''mu'' should be a single row vector.'); end
   if ~((size(mu, 2) == 1) || (size(mu, 2) == m)), error('Argument ''mu'' should be either a scalar or an m-elements row vector.'); end
   % SIGMA
   if isempty(sigma), error('Argument ''sigma'' should not be empty.'); end
   if ~isnumeric(sigma), error('Argument ''sigma'' should be a numeric value.'); end
   if (size(sigma, 1) > 1), error('Argument ''sigma'' should be a single row vector.'); end
   if ~((size(sigma, 2) == 1) || (size(sigma, 2) == m)), error('Argument ''sigma'' should be either a scalar or an m-elements row vector.'); end
   if any(sigma <= 0), error('Argument ''sigma'' should contain only positive values.'); end
   % L
   if isempty(l), error('Argument ''l'' should not be empty.'); end
   if ~isnumeric(l), error('Argument ''l'' should be a numeric value.'); end
   if (size(l, 1) > 1), error('Argument ''l'' should be a single row vector.'); end
   if ~((size(l, 2) == 1) || (size(l, 2) == m)), error('Argument ''l'' should be either a scalar or an m-elements row vector.'); end
   % U
   if isempty(u), error('Argument ''u'' should not be empty.'); end
   if ~isnumeric(u), error('Argument''u'' should be a numeric value.'); end
   if (size(u, 1) > 1), error('Argument ''u'' should be a single row vector.'); end
   if ~((size(u, 2) == 1) || (size(u, 2) == m)), error('Argument ''u'' should be either a scalar or an m-elements row vector.'); end
   if any(u <= l), error('All values of argument ''u'' must be higher than the corresponding values of ''l'''); end
   
   % Computing the boundaries probabilities
   p_l = normcdf((l - mu)./sigma);
   p_u = normcdf((u - mu)./sigma);
   % Computing
   R = rand(n, m);
   P = ones(n, m)*diag(p_l) + R*(diag(p_u) - diag(p_l));
   Z = ones(n, m)*diag(mu) + norminv(P)*diag(sigma);
end


