function val = wahi_ray_cache(key, builder)
% wahi_ray_cache  Memoise expensive ray-model objects (fans, eigenray paths,
% first-arrival tables) by a string key.
%
%   val = wahi_ray_cache(key, @() buildSomething())
%   wahi_ray_cache('clear')
%
% Each MATLAB worker (parfor) keeps its own cache. The cache is cleared when
% it exceeds 2000 entries.
persistent M
if isempty(M), M = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
if nargin == 1 && strcmp(key, 'clear')
    M = containers.Map('KeyType', 'char', 'ValueType', 'any');
    val = [];
    return;
end
if isKey(M, key)
    val = M(key);
    return;
end
val = builder();
if M.Count > 2000
    M = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
M(key) = val;
end
