% Kernel lineal para gp_reg.m
% where a is of size D by n, b is of size D by m (or empty), C and Q are of
% size n by m and c is of size D by 1.
%
% TFG Aero Rocío Navarro Villarino

function C = kernel_lin(a, b);

if nargin < 1 | nargin > 3 | nargout > 1
  error('Wrong number of arguments.');
end

if nargin == 1 | isempty(b)                   % input arguments are taken to be
  b = a;                                   % identical if b is missing or empty
end 

[D, n] = size(a); 
[d, m] = size(b);
if d ~= D
  error('Error: column lengths must agree.');
end

 C = zeros(n,m);
 for d = 1:D
%   C = C + (repmat(b(d,:), n, 1) - param).*(repmat(a(d,:)', 1, m) - param);
% %   C = C + repmat(b(d,:), n, 1).*repmat(a(d,:)', 1, m) - param;
  C = C + repmat(b(d,:), n, 1).*repmat(a(d,:)', 1, m);
 end
end

