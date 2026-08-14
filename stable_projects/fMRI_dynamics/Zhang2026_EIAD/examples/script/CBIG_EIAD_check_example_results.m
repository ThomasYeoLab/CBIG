function CBIG_EIAD_check_example_results(output_struct, ref_file, tol)
% CBIG_EIAD_CHECK_EXAMPLE_RESULTS  Compare example results against the reference.
%
% Compares the results produced by CBIG_EIAD_example against the reference
% output field by field, and throws an error if any field differs by more than
% a numeric tolerance.
%
% Inputs:
%   output_struct - struct of results produced by CBIG_EIAD_example.
%   ref_file      - (optional) path to the reference .mat file, which must hold
%                   a variable named 'output_struct'. Default:
%                   '../reference_output/reference_output.mat'.
%   tol           - (optional) absolute tolerance for the comparison.
%                   Default: 1e-6.
%
% On success, prints:
%   Passed! Results are the same as the reference outputs.
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 2 || isempty(ref_file)
    ref_file = '../reference_output/reference_output.mat';
end
if nargin < 3 || isempty(tol)
    tol = 1e-6;
end

reference_output = load(ref_file);
reference_output = reference_output.output_struct;

f_user = fieldnames(output_struct);

for i = 1:numel(f_user)
    x = output_struct.(f_user{i});
    y = reference_output.(f_user{i});

    % the GLM field is a table; compare it as a numeric array
    if istable(x), x = table2array(x); end
    if istable(y), y = table2array(y); end

    if ~isnumeric(x) || ~isnumeric(y) || ~isequal(size(x), size(y)) || ...
            any(abs(x(:) - y(:)) > tol)
        error(['Failed. Field ''%s'' differs from the reference output by more ' ...
               'than the tolerance of %g. Please check the hardware and ' ...
               'software versions.'], f_user{i}, tol)
    end
end

disp('Passed! Results are the same as the reference outputs.')

end
