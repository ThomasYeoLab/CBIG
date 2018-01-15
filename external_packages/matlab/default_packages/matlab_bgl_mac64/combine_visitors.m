function cv = combine_visitors(varargin)
% COMBINE_VISITORS Generate a new visitor by combining existing visitors
%
% cv = combine_visitors(v1, v2, ...) generates a new algorithm visitor that
% works by calling the functions in v1 followed by the functions in v2, and
% so on.  
%
% The value returned by the combined visitor function is the bitwise & of
% all the individual return values.  So, if any visitor requests the
% algorithm to halt, then the algorithm will halt.
%
% Note: using a combined visitor is somewhat slower than writing a custom
% combined visitor yourself.
%
% Example:
%    vis1 = struct();
%    vis1.examine_vertex = @(u) fprintf('vis1: examine_vertex(%i)\n', u);
%    vis2 = struct();
%    vis2.examine_vertex = @(u) fprintf('vis2: examine_vertex(%i)\n', u);
%    combined_vis = combine_visitors(vis1, vis2);
%    load graphs/bfs_example.mat
%    breadth_first_search(A,1,combined_vis);

% David Gleich
% Copyright, Stanford University, 2006-2008

% trivial output
if isempty(varargin)
    error('matlab_bgl:invalidParameter', 'combine_visitors requires at least one argument');
end

if length(varargin) == 1, cv = varargin{1}; return; end

cv_fn = struct();

for ii = 1:length(varargin)
    fn = fieldnames(varargin{ii});
    
    for jj = 1:length(fn)
        if ~isfield(cv_fn,fn{jj})
            cv_fn.(fn{jj}) = 1;
        else
            cv_fn.(fn{jj}) = cv_fn.(fn{jj}) + 1;
        end
    end
end

    function cfunc=combine_function(name,count)
        
        used_visitors = cell(count,1);
        cur_used = 1;
        
        for kk = 1:length(varargin)
            if isfield(varargin{kk}, name)
                used_visitors{cur_used} = varargin{kk}.(name);
                cur_used = cur_used + 1;
            end
        end
        
        function rval = combine_function_impl(varargin)
            rval = 1;
            for ll=1:length(used_visitors)
                try
                    stop = feval(used_visitors{ll}, varargin{:});
                    rval = rval & double(stop);
                catch
                    [lastmsg,lastid] = lasterr;
                    if ~strcmp(lastid, 'MATLAB:TooManyOutputs')
                        rethrow(lasterr);
                    else
                        feval(used_visitors{ll}, varargin{:});
                    end
                end
            end
        end
        
        cfunc = @combine_function_impl;
        
    end


cv = struct();

% specify only the visitors functions that only occured once
for ii = 1:length(varargin)
    fn = fieldnames(varargin{ii});
    
    for jj = 1:length(fn)
        if (cv_fn.(fn{jj}) == 1)
            cv.(fn{jj}) = varargin{ii}.(fn{jj});
        end
    end
end

% specify all visitors that occured more than once
fn = fieldnames(cv_fn);
for jj=1:length(fn)
    if (cv_fn.(fn{jj}) > 1)
        cv.(fn{jj}) = combine_function(fn{jj}, cv_fn.(fn{jj}));
    end
end

% overall function end
end
