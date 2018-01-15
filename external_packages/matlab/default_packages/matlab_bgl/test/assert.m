function assert(condition,varargin)

if condition, return;
else
  error(varargin{:});
end
