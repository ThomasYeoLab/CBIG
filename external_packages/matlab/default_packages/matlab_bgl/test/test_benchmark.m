function out=test_benchmark

%
% David Gleich
% Copyright, Stanford University, 2007
%

% 
% 9 July 2007
% Initial version
%

%% History

%% 
%  2008-10-07, Version 2.1, Matlab 2007b, boost 1.33.0, 
%    g++-3.4 (lib), gcc-? (mex)
%         airfoil       west    cs-stan    minneso      tapir   
%  large   0.223 s    0.024 s    0.390 s    0.073 s    0.046 s  
%    med     NaN s    0.955 s      NaN s      NaN s    6.621 s  
%  small     NaN s    0.758 s      NaN s      NaN s      NaN s  

%% 
%  2008-10-07: Version 3.0, Matlab 2007b, boost 1.34.1, 
%    g++-4.0 (lib), gcc-? (mex)
%
%         airfoil       west    cs-stan    minneso      tapir   
%  large   0.183 s    0.017 s    0.222 s    0.048 s    0.037 s  
%    med     NaN s    0.593 s      NaN s      NaN s    3.901 s  
%  small     NaN s    0.543 s      NaN s      NaN s      NaN s  

%% 
%  2008-10-07: Version 4.0, Matlab 2007b, Boost 1.36.0,
%    g++-3.4 (lib) and gcc-4.1 (mex)
%
%         airfoil       west    cs-stan    minneso      tapir   
%  large   0.151 s    0.019 s    0.260 s    0.054 s    0.031 s  
%    med     NaN s    0.673 s      NaN s      NaN s    4.459 s  
%  small     NaN s    0.665 s      NaN s      NaN s      NaN s 

%% Benchmarks for MatlabBGL
% This function computes a series of benchmarks 

%% Create tests
% There are three types of tests, large scale, medium scale, and small
% scale.  
nrep = 10;
rand('state',0);

%% Setup results output
r.dt_large = 0;
r.dt_med = 0;
r.dt_small = 0;
testi = 1;

%% airfoil mesh
% Forming the sparse adjacency matrix and making it positive definite
af = load('airfoil');
n = max(max(af.i),max(af.j));
A = sparse(af.i,af.j,-1,n,n);
A = A + A';
d = abs(sum(A)) + 1;
A = A + diag(sparse(d));

r.name = 'airfoil'; fprintf('Graph: %s\n', r.name);
r.dt_large = bench_large(A,nrep);
r.dt_med = NaN;
r.dt_small = NaN;

results(testi) = r;
testi = testi+1;

%% west0479 matrix
load west0479;
A = west0479;

r.name = 'west'; fprintf('Graph: %s\n', r.name);
r.dt_large = bench_large(A,nrep);
r.dt_med = bench_medium(A,nrep);
r.dt_small = bench_small(A,nrep);

results(testi) = r;
testi = testi+1;

%% cs-stanford matrix
load('../graphs/cs-stanford');

r.name = 'cs-stan'; fprintf('Graph: %s\n', r.name);
r.dt_large = bench_large(A,nrep);
r.dt_med = NaN;
r.dt_small = NaN;

results(testi) = r;
testi = testi+1;

%% minnesota
load('../graphs/minnesota');

r.name = 'minneso'; fprintf('Graph: %s\n', r.name);
r.dt_large = bench_large(A,nrep);
r.dt_med = NaN;
r.dt_small = NaN;

results(testi) = r;
testi = testi+1;

%% tapir
load('../graphs/tapir');

r.name = 'tapir'; fprintf('Graph: %s\n', r.name);
r.dt_large = bench_large(A,nrep);
r.dt_med = bench_medium(A,nrep);
r.dt_small = NaN;

results(testi) = r;
testi = testi+1;

%% display results
fprintf('%5s ', ''); cellfun(@(s) fprintf(' %7s   ', s), {results.name}); fprintf('\n');
fprintf('%5s ','large'); fprintf('%7.3f s  ', [results.dt_large]); fprintf('\n');
fprintf('%5s ','med');   fprintf('%7.3f s  ', [results.dt_med]); fprintf('\n');
fprintf('%5s ','small');   fprintf('%7.3f s  ', [results.dt_small]); fprintf('\n');

if nargout > 1
    out = results;
end

end

    function dt=bench_large(A,nrep)
        As = spones(A);
        Asym = A+A';

        dt = 0;
        tic;
        for ri=1:nrep
            c=components(A);
        end
        op_dt = toc;
        dt = dt + op_dt;

        fprintf('  %30s: %7.3f s\n', 'components', op_dt);

        dt = 0;
        tic;
        for ri=1:nrep
            c=biconnected_components(Asym);
        end
        op_dt = toc;
        dt = dt + op_dt;

        fprintf('  %30s: %7.3f s\n', 'biconnect_components', op_dt);
        
        tic;
        for ri=1:nrep
            u = ceil(rand(1)*size(A,1)); v = ceil(rand(1)*size(A,1));
            f=max_flow(As,u,v);
        end
        op_dt = toc;
        dt = dt + op_dt;

        fprintf('  %30s: %7.3f s\n', 'max_flow', op_dt);
 
        tic;
        for ri=1:nrep
            u = ceil(rand(1)*size(A,1)); v = ceil(rand(1)*size(A,1));
            [i j v] = mst(Asym);
        end
        op_dt = toc;
        dt = dt + op_dt;

        fprintf('  %30s: %7.3f s\n', 'mst', op_dt);
    end

    function dt=bench_medium(A,nrep)
        dt = 0;
        tic;
        for ri=1:nrep
            bc=betweenness_centrality(spones(A));
        end
        op_dt = toc;
        dt = dt + op_dt;

        fprintf('  %30s: %7.3f s\n', 'betweenness_centrality', op_dt);
    end

    function dt=bench_small(A,nrep)
        dt = 0;
        tic;
        for ri=1:nrep
            bc=all_shortest_paths(spones(A));
        end
        op_dt = toc;
        dt = dt + op_dt;

        fprintf('  %30s: %7.3f s\n', 'all_shortest_paths', op_dt);
    end
