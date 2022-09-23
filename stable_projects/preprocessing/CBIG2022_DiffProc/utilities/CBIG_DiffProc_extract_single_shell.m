function CBIG_extract_single_shell(outdir, subjID, bval_file,bvec_file, shell, round_1000)
    % Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    % extract single shell and b0 image from bvals and bvec in fsl format
    bvals = dlmread(bval_file);
    if round_1000 == 'y'
        bvals = round(bvals/1000)*1000;
    end
    bvecs = dlmread(bvec_file);

    new_idx = 0;
    for n = 1:length(bvals)
        if (bvals(n) == 0 || bvals(n) == shell)
            new_idx = new_idx + 1;
            %frame_loc(new_idx) = n;
            new_bvals(new_idx) = bvals(n);
            new_bvecs(:,new_idx) = bvecs(:,n);
        end
    end
       
    % save modified bvals and bvecs
    dlmwrite(fullfile(outdir, strcat(subjID, '_single_shell_', num2str(shell) ,'.bval')), new_bvals,'delimiter', ' ');
    dlmwrite(fullfile(outdir, strcat(subjID, '_single_shell_', num2str(shell) ,'.bvec')), new_bvecs,'delimiter', ' ');
end