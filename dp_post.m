% Input:
%     data         a row vector. The sampled positions according to some
%                  distribution.
%     alpha        a scalar. The concentration parameter of the DP.
%     actN         a scalar. The number of the maximum of activated atoms.
%     maxIter      a scalar. The maximum number of gibbs iterations.
% Output:
%     distro_mean  a row vector of length actN.
function distro_mean = dp_post(data, alpha, actN, maxIter)
if nargin < 3
    actN = 100;
end
if nargin < 4
    maxIter = 500;
end

%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------
distro_mean = zeros(1, actN);

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------
for iter = 1:maxIter
    % a_k = m_k
    a = histcounts(data, 1:actN+1);
    % b_k = m_{k+1} + ... + m_{infty}
    b = cumsum(a, 'reverse');
    b = b(2:end);
    b = [b, 0];
    
    % prior is beta(1, alpha)
    a = a + 1;
    b = b + alpha;
    
    % sample V
    V = betarnd(a, b);
    
    % compute distro
    distro = V;
    V = cumprod(1 - V);
    distro(2:end) = distro(2:end) .* V(1:end-1);
    
    % acculate
    distro_mean = distro + distro_mean;
    
end
distro_mean = distro_mean / maxIter;
end