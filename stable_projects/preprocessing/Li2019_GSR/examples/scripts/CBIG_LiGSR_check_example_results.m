function CBIG_LiGSR_check_example_results

% CBIG_LiGSR_check_example_results
% 
% This function compares the example tests results with the reference
% results.
%
% Inputs: this function has no input variable
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Jingwei Li

fprintf('Comparing example results ...\n')
exm_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
    'Li2019_GSR', 'examples');
outdir = fullfile(exm_dir, 'output');
refdir = fullfile(exm_dir, 'ref_output');
nbehaviors = 2;

%% compare results of variance component model
jk_seeds = 10;
d = 2;
fprintf('1. Variance Component model:\n')
for seed = 1:jk_seeds
    fprintf('\tJackknife #%d ...\n', seed)
    for b = 1:nbehaviors
        fprintf('\t\tBehavior %d ...\n', b)
        
        outfile = fullfile(outdir, 'VarianceComponentModel', ['del' num2str(d) '_set' ...
            num2str(seed)], ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
        reffile = fullfile(refdir, 'VarianceComponentModel', ['del' num2str(d) '_set' ...
            num2str(seed)], ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
        out = load(outfile);
        ref = load(reffile);
        
        reffields = fieldnames(ref.morpho);
        outfields = fieldnames(out.morpho);
        assert(length(reffields) == length(outfields), ...
            '\t\t\tOutput ''morpho'' has wrong structure.')
        
        for i = 1:length(reffields)
            curr_outfield = getfield(out.morpho, outfields{i});
            curr_reffield = getfield(ref.morpho, reffields{i});
            assert(isequal(size(curr_reffield), size(curr_outfield)), ...
                sprintf('\t\t\tstructure field %s is of wrong size.', reffields{i}))
            
            if(all(isnan(curr_reffield)))
                assert(all(isnan(curr_outfield)), sprintf(...
                    '\t\t\treference field %s is NaN, but output field %s is not NaN.\n', ...
                    reffields{i}, reffields{i}))
            else
                maxdif = max(abs(curr_reffield(:) - curr_outfield(:)));
                assert(maxdif<1e-10, sprintf(...
                    '\t\t\tstructure field %s differed by (max abs diff) %f.', ...
                    reffields{i}, maxdif));
            end
        end
    end
end

fprintf('\tFull set ...\n')
for b = 1:nbehaviors
    fprintf('\t\tBehavior %d ...\n', b)
    
    outfile = fullfile(outdir, 'VarianceComponentModel', 'fullset', ...
        ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
    reffile = fullfile(refdir, 'VarianceComponentModel', 'fullset', ...
        ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
    out = load(outfile);
    ref = load(reffile);
    
    reffields = fieldnames(ref.morpho);
    outfields = fieldnames(out.morpho);
    assert(length(reffields) == length(outfields), '\t\t\tOutput ''morpho'' has wrong structure.')
    
    for i = 1:length(reffields)
        curr_outfield = getfield(out.morpho, outfields{i});
        curr_reffield = getfield(ref.morpho, reffields{i});
        assert(isequal(size(curr_reffield), size(curr_outfield)), ...
            sprintf('\t\t\tstructure field %s is of wrong size.', reffields{i}))
        
        if(all(isnan(curr_reffield)))
            assert(all(isnan(curr_outfield)), sprintf(...
                '\t\t\trefernce field %s is NaN, but output field %s is not NaN.\n', ...
                reffields{i}, reffields{i}))
        else
            maxdif = max(abs(curr_reffield(:) - curr_outfield(:)));
            assert(maxdif<1e-8, sprintf(...
                '\t\t\tstructure field %s differed by (max abs diff) %f.', ...
                reffields{i}, maxdif));
        end
    end
end
fprintf('Variance component model results were replicated.\n')


%% compare results of kernel regression
fprintf('2. Kernel regression:\n')
outfile = fullfile(outdir, 'KernelRidgeRegression', 'final_result.mat');
reffile = fullfile(refdir, 'KernelRidgeRegression', 'final_result.mat');
out = load(outfile);
ref = load(reffile);

out = orderfields(out,ref);
reffields = fieldnames(ref);
outfields = fieldnames(out);
assert(length(reffields) == length(outfields), 'Output has wrong variables.')

for i = 1:length(reffields)
    curr_outfield = getfield(out, outfields{i});
    curr_reffield = getfield(ref, reffields{i});
    assert(isequal(size(curr_reffield), size(curr_outfield)), ...
        sprintf('\t\t\tvariable %s is of wrong size.', reffields{i}));
    
    if(strcmp(reffields{i}, 'optimal_stats') || strcmp(reffields{i}, 'optimal_kernel'))
        subfields = fieldnames(curr_reffield);
        for j = 1:length(subfields)
            sub_reffield = getfield(curr_reffield, subfields{j});
            sub_outfield = getfield(curr_outfield, subfields{j});
            
            if(all(reshape(isnan(sub_reffield), [numel(sub_reffield) 1])))
                assert(all(reshape(isnan(sub_outfield), [numel(sub_outfield) 1])), sprintf(...
                    '\treference variable %s.%s is Nan, but output is not.\n', ...
                    reffields{i}, subfields{j}));
            else
                maxdif = max(abs(sub_reffield(:) - sub_outfield(:)));
                assert(maxdif<1e-10, sprintf('\tvariable %s.%s differed by (max abs diff) %f', ...
                    reffields{i}, subfields{j}, maxdif))
            end
        end
    else
        if(all(reshape(isnan(curr_reffield), [numel(curr_reffield) 1])))
            assert(all(reshape(isnan(curr_outfield), [numel(curr_outfield) 1])), sprintf(...
                '\treference variable %s is NaN, but output is not.\n', reffields{i}))
        else
            maxdif = max(abs(curr_reffield(:) - curr_outfield(:)));
            assert(maxdif<1e-10, sprintf('\tvariable %s differed by (max abs diff) %f.', ...
                reffields{i}, maxdif));
        end
    end
end
fprintf('Kernel regression results were replicated.\n')

%% remove output
rmdir(outdir, 's')

end

