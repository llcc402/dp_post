% Input:
%     distro    a row vector. It is generated by a Dirichlet Process.
%     n         a scalar. The size of the data set.
% Output:
%     data      a row vector of length n. 
function data = data_generate(distro, n)
if nargin < 1
    distro = gem(100, 5);
end
if nargin < 2
    n = 500;
end

r = rand(1, n);
d = [0, cumsum(distro)];
[~, ~, data] = histcounts(r, d);
end