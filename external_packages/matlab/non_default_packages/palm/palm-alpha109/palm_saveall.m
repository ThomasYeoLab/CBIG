function palm_saveall(plm,opts)
% Save most of the outputs from PALM.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2014
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% For the MV and NPC methods in which the most significant stats are the
% smallest, rather than the largest, use reverse comparisons.
if opts.NPC,
    if plm.npcrev,
        npcextr = @min;
    else
        npcextr = @max;
    end
end

% For the "noperm" approximation, the statistic didn't have to be saved
% that early, since it's all pretty fast. Save it now then.
if opts.accel.noperm,
    if opts.saveunivariate,
        for y = 1:plm.nY,
            if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
            for m = loopM,
                for c = 1:plm.nC(m),
                    palm_quicksave(plm.G{y}{m}{c},0,opts,plm,y,m,c, ...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                end
            end
        end
    end
    if opts.MV,
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave(plm.Q{m}{c},0,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{m}{c}));
            end
        end
    end
end


% Start with the uncorrected, but don't save them yet.
% They'll be used later for the FDR.
fprintf('Computing p-values.\n');
if opts.saveuncorrected,
    for y = 1:plm.nY,
        if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
        for m = loopM,
            for c = 1:plm.nC(m),
                if opts.accel.negbin,
                    idx = plm.Gpperm{y}{m}{c} == opts.accel.negbin;
                    plm.Gpperm{y}{m}{c}( idx) = (plm.Gpperm{y}{m}{c}( idx)    )./(plm.Gppermp{y}{m}{c}( idx)); % Haldane (1945), but need to add 1 to num and den, see Besag & Clifford (1991).
                    plm.Gpperm{y}{m}{c}(~idx) = (plm.Gpperm{y}{m}{c}(~idx) + 1)./plm.nP{m}(c); % the unpermuted isn't counted here, hence the +1. Also, this is the same as if dividing by (plm.Gppermp{y}{m}{c}(~idx) + 1)
                elseif opts.accel.tail,
                    for t = 1:plm.Ysiz(y),
                        plm.Gpperm{y}{m}{c}(1,t) = palm_pareto(plm.G{y}{m}{c}(1,t),plm.Gperms{y}{m}{c}(:,t),false,opts.accel.tail_thr,opts.accel.G1out);
                    end
                elseif opts.accel.gamma,
                    for t = 1:plm.Ysiz(y),
                        plm.Gpperm{y}{m}{c}(1,t) = approxgamma(plm.G{y}{m}{c}(1,t),plm.Gperms{y}{m}{c}(:,t),false,1/plm.nP{m}(c),opts.accel.G1out);
                    end
                elseif ~ opts.accel.noperm, % the "noperm" case is already treated
                    plm.Gpperm{y}{m}{c} = plm.Gpperm{y}{m}{c}/plm.nP{m}(c);
                end
            end
        end
    end
    if opts.tfce.uni.do,
        for y = 1:plm.nY,
            if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
            for m = loopM,
                for c = 1:plm.nC(m),
                    if opts.accel.tail,
                        for t = 1:plm.Ysiz(y),
                            plm.Gtfcepperm{y}{m}{c}(1,t) = palm_pareto(plm.Gtfce{y}{m}{c}(1,t),plm.Gtfceperms{y}{m}{c}(:,t),false,opts.accel.tail_thr,opts.accel.G1out);
                        end
                    elseif opts.accel.gamma,
                        for t = 1:plm.Ysiz(y),
                            plm.Gtfcepperm{y}{m}{c}(1,t) = approxgamma(plm.Gtfce{y}{m}{c}(1,t),plm.Gtfceperms{y}{m}{c}(:,t),false,1/plm.nP{m}(c),opts.accel.G1out);
                        end
                    else
                        plm.Gtfcepperm{y}{m}{c} = plm.Gtfcepperm{y}{m}{c}/plm.nP{m}(c);
                    end
                end
            end
        end
    end
    if opts.NPC,
        if opts.npcmod && ~ opts.npccon,
            if opts.designperinput, loopM = 1; else loopM = 1:plm.nM; end
            for m = loopM,
                for c = 1:plm.nC(m),
                    if opts.accel.tail,
                        for t = 1:plm.Ysiz(1),
                            plm.Tpperm{m}{c}(1,t) = palm_pareto(plm.T{m}{c}(1,t),plm.Tperms{m}{c}(:,t),plm.npcrev,opts.accel.tail_thr,opts.accel.G1out);
                        end
                    elseif opts.accel.gamma,
                        for t = 1:plm.Ysiz(1),
                            plm.Tpperm{m}{c}(1,t) = approxgamma(plm.T{m}{c}(1,t),plm.Tperms{m}{c}(:,t),plm.npcrev,1/plm.nP{m}(c),opts.accel.G1out);
                        end
                    else
                        plm.Tpperm{m}{c} = plm.Tpperm{m}{c}/plm.nP{m}(c);
                    end
                end
            end
            if opts.tfce.npc.do,
                for m = loopM,
                    for c = 1:plm.nC(m),
                        if opts.accel.tail,
                            for t = 1:plm.Ysiz(1),
                                plm.Ttfcepperm{m}{c}(1,t) = palm_pareto(plm.Ttfce{m}{c}(1,t),plm.Ttfceperms{m}{c}(:,t),false,opts.accel.tail_thr,opts.accel.G1out);
                            end
                        elseif opts.accel.gamma,
                            for t = 1:plm.Ysiz(1),
                                plm.Ttfcepperm{m}{c}(1,t) = approxgamma(plm.Ttfce{m}{c}(1,t),plm.Ttfceperms{m}{c}(:,t),false,1/plm.nP{m}(c),opts.accel.G1out);
                            end
                        else
                            plm.Ttfcepperm{m}{c} = plm.Ttfcepperm{m}{c}/plm.nP{m}(c);
                        end
                    end
                end
            end
        elseif opts.npccon,
            for j = 1:numel(plm.Tmax),
                if opts.accel.tail,
                    for t = 1:plm.Ysiz(1),
                        plm.Tpperm{j}(1,t) = palm_pareto(plm.T{j}(1,t),plm.Tperms{j}(:,t),plm.npcrev,opts.accel.tail_thr,opts.accel.G1out);
                    end
                elseif opts.accel.gamma,
                    for t = 1:plm.Ysiz(1),
                        plm.Tpperm{j}(1,t) = approxgamma(plm.T{j}(1,t),plm.Tperms{j}(:,t),plm.npcrev,1/plm.nP{1}(1),opts.accel.G1out);
                    end
                else
                    plm.Tpperm{j} = plm.Tpperm{j}/size(plm.Tmax{j},1);
                end
            end
            if opts.tfce.npc.do,
                for j = 1:numel(plm.Tmax),
                    if opts.accel.tail,
                        for t = 1:plm.Ysiz(1),
                            plm.Ttfcepperm{j}(1,t) = palm_pareto(plm.Ttfce{j}(1,t),plm.Ttfceperms{j}(:,t),false,opts.accel.tail_thr,opts.accel.G1out);
                        end
                    elseif opts.accel.gamma,
                        for t = 1:plm.Ysiz(1),
                            plm.Ttfcepperm{j}(1,t) = approxgamma(plm.Ttfce{j}(1,t),plm.Ttfceperms{j}(:,t),false,1/plm.nP{1}(1),opts.accel.G1out);
                        end
                    else
                        plm.Ttfcepperm{j} = plm.Ttfcepperm{j}/size(plm.Ttfcemax{j},1);
                    end
                end
            end
        end
    end
    if opts.MV || opts.CCA || opts.PLS,
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                if opts.accel.negbin,
                    idx = plm.Qpperm{m}{c} == opts.accel.negbin;
                    plm.Qpperm{m}{c}( idx) = (plm.Qpperm{m}{c}( idx) - 1)./(plm.Qppermp{m}{c}( idx) - 1); % Haldane (1945)
                    plm.Qpperm{m}{c}(~idx) = (plm.Qpperm{m}{c}(~idx) + 1)./plm.nP{m}(c); % same as if dividing by (plm.Qppermp{m}{c}(~idx) + 1)
                elseif opts.accel.tail,
                    for t = 1:plm.Ysiz(1),
                        plm.Qpperm{m}{c}(1,t) = palm_pareto(plm.Q{m}{c}(1,t),plm.Qperms{m}{c}(:,t),plm.mvrev{m}{c},opts.accel.tail_thr,opts.accel.G1out);
                    end
                elseif opts.accel.gamma,
                    for t = 1:plm.Ysiz(1),
                        plm.Qpperm{m}{c}(1,t) = approxgamma(plm.Q{m}{c}(1,t),plm.Qperms{m}{c}(:,t),plm.mvrev{m}{c},1/plm.nP{m}(c),opts.accel.G1out);
                    end
                elseif ~ opts.accel.noperm, % the "noperm" case is already treated
                    plm.Qpperm{m}{c} = plm.Qpperm{m}{c}/plm.nP{m}(c);
                end
            end
        end
        if opts.tfce.mv.do,
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    if opts.accel.tail,
                        for t = 1:plm.Ysiz(1),
                            plm.Qtfcepperm{m}{c}(1,t) = palm_pareto(plm.Qtfce{m}{c}(1,t),plm.Qtfceperms{m}{c}(:,t),false,opts.accel.tail_thr,opts.accel.G1out);
                        end
                    elseif opts.accel.gamma,
                        for t = 1:plm.Ysiz(1),
                            plm.Qtfcepperm{m}{c}(1,t) = approxgamma(plm.Qtfce{m}{c}(1,t),plm.Qtfceperms{m}{c}(:,t),false,1/plm.nP{m}(c),opts.accel.G1out);
                        end
                    else
                        plm.Qtfcepperm{m}{c} = plm.Qtfcepperm{m}{c}/plm.nP{m}(c);
                    end
                end
            end
        end
    end
end

% Save uncorrected & FWER-corrected within modality for this contrast.
if opts.saveunivariate,
    fprintf('Saving p-values (uncorrected, and corrected within modality and within contrast).\n');
    for y = 1:plm.nY,
        if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
        for m = loopM,
            for c = 1:plm.nC(m),
                
                % Only permutation p-value and its FDR ajustment are saved in the negative binomial mode.
                if opts.accel.negbin,
                    
                    % Permutation p-value, uncorrected (the -nouncorrected
                    % is caught in palm_takeargs.m already)
                    palm_quicksave(plm.Gpperm{y}{m}{c},1,opts,plm,y,m,c, ...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_uncp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    
                    % Permutation p-value, FDR adjusted
                    if opts.FDR,
                        palm_quicksave(fastfdr(plm.Gpperm{y}{m}{c}),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                else
                    
                    % Permutation p-value
                    if opts.saveuncorrected,
                        palm_quicksave(plm.Gpperm{y}{m}{c},1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_uncp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                    
                    % FWER-corrected
                    if opts.accel.tail,
                        Ptosave = palm_pareto(plm.G{y}{m}{c},plm.Gmax{y}{m}{c},false,opts.accel.tail_thr,opts.accel.G1out);
                    elseif opts.accel.noperm,
                        Ptosave = [];
                    elseif opts.accel.gamma,
                        Ptosave = approxgamma(plm.G{y}{m}{c},plm.Gmax{y}{m}{c},false,1/plm.nP{m}(c),opts.accel.G1out);
                    else
                        Ptosave = palm_datapval(plm.G{y}{m}{c},plm.Gmax{y}{m}{c},false);
                    end
                    palm_quicksave( ...
                        Ptosave,1,opts,plm,y,m,c,...
                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    
                    % Permutation p-value, FDR adjusted
                    if opts.FDR,
                        palm_quicksave(fastfdr(plm.Gpperm{y}{m}{c}),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                    
                    % Cluster statistic results.
                    if opts.cluster.uni.do,
                        
                        % Cluster statistic.
                        palm_quicksave(plm.Gclu{y}{m}{c},0,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,opts.cluster.str,plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        
                        % Cluster statistic FWER p-value
                        if opts.accel.tail,
                            Ptosave = palm_pareto(plm.Gclu{y}{m}{c},plm.Gclumax{y}{m}{c},false,opts.accel.tail_thr,opts.accel.G1out);
                        elseif opts.accel.gamma,
                            Ptosave = approxgamma(plm.Gclu{y}{m}{c},plm.Gclumax{y}{m}{c},false,1/plm.nP{m}(c),opts.accel.G1out);
                        else
                            Ptosave = palm_datapval(plm.Gclu{y}{m}{c},plm.Gclumax{y}{m}{c},false);
                        end
                        palm_quicksave( ...
                            Ptosave,1,opts,plm,y,m,c,...
                            sprintf('%s',opts.o,opts.cluster.str,plm.Gname{m}{c},'_fwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                    
                    % TFCE results
                    if opts.tfce.uni.do,
                        
                        % TFCE statistic
                        palm_quicksave(plm.Gtfce{y}{m}{c},0,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        
                        % TFCE p-value
                        if opts.saveuncorrected,
                            palm_quicksave(plm.Gtfcepperm{y}{m}{c},1,opts,plm,y,m,c,...
                                sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_uncp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                        
                        % TFCE FWER-corrected within modality and contrast.
                        if opts.accel.tail,
                            Ptosave = palm_pareto(plm.Gtfce{y}{m}{c},plm.Gtfcemax{y}{m}{c},false,opts.accel.tail_thr,opts.accel.G1out);
                        elseif opts.accel.gamma,
                            Ptosave = approxgamma(plm.Gtfce{y}{m}{c},plm.Gtfcemax{y}{m}{c},false,1/plm.nP{m}(c),opts.accel.G1out);
                        else
                            Ptosave = palm_datapval(plm.Gtfce{y}{m}{c},plm.Gtfcemax{y}{m}{c},false);
                        end
                        palm_quicksave( ...
                            Ptosave,1,opts,plm,y,m,c,...
                            sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_fwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        
                        % TFCE p-value, FDR adjusted.
                        if opts.FDR,
                            palm_quicksave(fastfdr(plm.Gtfcepperm{y}{m}{c}),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_fdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                end
                
                % Parametric p-value and its FDR adjustment
                if opts.savepara,
                    if opts.saveuncorrected,
                        P = palm_quicksave(plm.G{y}{m}{c},2,opts,plm,y,m,c,...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_uncparap',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                    if opts.FDR,
                        palm_quicksave(fastfdr(P),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_fdrparap',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                end
            end
        end
    end
    
    % Save FWER & FDR corrected across modalities.
    if opts.corrmod,
        fprintf('Saving p-values (corrected across modalities).\n')
        
        % FDR correction (non-spatial stats)
        if opts.FDR,
            if opts.designperinput,
                for c = 1:plm.nC(1),
                    pmerged = zeros(sum(plm.Ysiz),1);
                    for y = 1:plm.nY,
                        m = y;
                        pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gpperm{y}{m}{c};
                    end
                    pfdradj = fastfdr(pmerged);
                    for y = 1:plm.nY,
                        m = y;
                        palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                end
            else
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        pmerged = zeros(sum(plm.Ysiz),1);
                        for y = 1:plm.nY,
                            pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gpperm{y}{m}{c};
                        end
                        pfdradj = fastfdr(pmerged);
                        for y = 1:plm.nY,
                            palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                end
            end
        end
        
        if ~ opts.accel.negbin,
            
            % FWER correction (non-spatial stats)
            if opts.designperinput,
                for c = 1:plm.nC(1),
                    distmax = zeros(plm.nP{1}(c),plm.nY);
                    for y = 1:plm.nY,
                        m = y;
                        distmax(:,y) = plm.Gmax{y}{m}{c};
                    end
                    distmax = max(distmax,[],2);
                    for y = 1:plm.nY,
                        m = y;
                        if opts.accel.tail,
                            Ptosave = palm_pareto(plm.G{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                        elseif opts.accel.gamma,
                            Ptosave = approxgamma(plm.G{y}{m}{c},distmax,false,1/plm.nP{1}(c),opts.accel.G1out);
                        else
                            Ptosave = palm_datapval(plm.G{y}{m}{c},distmax,false);
                        end
                        palm_quicksave( ...
                            Ptosave,1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                end
            else
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        distmax = zeros(plm.nP{m}(c),plm.nY);
                        for y = 1:plm.nY,
                            distmax(:,y) = plm.Gmax{y}{m}{c};
                        end
                        distmax = max(distmax,[],2);
                        for y = 1:plm.nY,
                            if opts.accel.tail,
                                Ptosave = palm_pareto(plm.G{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                            elseif opts.accel.gamma,
                                Ptosave = approxgamma(plm.G{y}{m}{c},distmax,false,1/plm.nP{m}(c),opts.accel.G1out);
                            else
                                Ptosave = palm_datapval(plm.G{y}{m}{c},distmax,false);
                            end
                            palm_quicksave( ...
                                Ptosave,1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                end
            end
            
            % Cluster statistic
            if opts.cluster.uni.do && ...
                    (all(plm.Yisvol) || all(plm.Yissrf)),
                if opts.designperinput,
                    for c = 1:plm.nC(1),
                        distmax = zeros(plm.nP{1}(c),plm.nY);
                        for y = 1:plm.nY,
                            m = y;
                            distmax(:,y) = plm.Gclumax{y}{m}{c};
                        end
                        distmax = max(distmax,[],2);
                        for y = 1:plm.nY,
                            m = y;
                            if opts.accel.tail,
                                Ptosave = palm_pareto(plm.Gclu{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                            elseif opts.accel.gamma,
                                Ptosave = approxgamma(plm.Gclu{y}{m}{c},distmax,false,1/plm.nP{1}(c),opts.accel.G1out);
                            else
                                Ptosave = palm_datapval(plm.Gclu{y}{m}{c},distmax,false);
                            end
                            palm_quicksave( ...
                                Ptosave,1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,opts.cluster.str,plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                else
                    for m = 1:plm.nM,
                        for c = 1:plm.nC(m),
                            distmax = zeros(plm.nP{m}(c),plm.nY);
                            for y = 1:plm.nY,
                                distmax(:,y) = plm.Gclumax{y}{m}{c};
                            end
                            distmax = max(distmax,[],2);
                            for y = 1:plm.nY,
                                if opts.accel.tail,
                                    Ptosave = palm_pareto(plm.Gclu{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                                elseif opts.accel.gamma,
                                    Ptosave = approxgamma(plm.Gclu{y}{m}{c},distmax,false,1/plm.nP{m}(c),opts.accel.G1out);
                                else
                                    Ptosave = palm_datapval(plm.Gclu{y}{m}{c},distmax,false);
                                end
                                palm_quicksave( ...
                                    Ptosave,1,opts,plm,y,m,c, ...
                                    sprintf('%s',opts.o,opts.cluster.str,plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                            end
                        end
                    end
                end
            end
            
            % TFCE
            if opts.tfce.uni.do && ...
                    (all(plm.Yisvol) || all(plm.Yissrf)),
                if opts.designperinput,
                    for c = 1:plm.nC(1),
                        distmax = zeros(plm.nP{1}(c),plm.nY);
                        for y = 1:plm.nY,
                            m = y;
                            distmax(:,y) = plm.Gtfcemax{y}{m}{c};
                        end
                        distmax = max(distmax,[],2);
                        for y = 1:plm.nY,
                            m = y;
                            if opts.accel.tail,
                                Ptosave = palm_pareto(plm.Gtfce{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                            elseif opts.accel.gamma,
                                Ptosave = approxgamma(plm.Gtfce{y}{m}{c},distmax,false,1/plm.nP{1}(c),opts.accel.G1out);
                            else
                                Ptosave = palm_datapval(plm.Gtfce{y}{m}{c},distmax,false);
                            end
                            palm_quicksave( ...
                                Ptosave,1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                else
                    for m = 1:plm.nM,
                        for c = 1:plm.nC(m),
                            distmax = zeros(plm.nP{m}(c),plm.nY);
                            for y = 1:plm.nY,
                                distmax(:,y) = plm.Gtfcemax{y}{m}{c};
                            end
                            distmax = max(distmax,[],2);
                            for y = 1:plm.nY,
                                if opts.accel.tail,
                                    Ptosave = palm_pareto(plm.Gtfce{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                                elseif opts.accel.gamma,
                                    Ptosave = approxgamma(plm.Gtfce{y}{m}{c},distmax,false,1/plm.nP{m}(c),opts.accel.G1out);
                                else
                                    Ptosave = palm_datapval(plm.Gtfce{y}{m}{c},distmax,false);
                                end
                                palm_quicksave( ...
                                    Ptosave,1,opts,plm,y,m,c, ...
                                    sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_mfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                            end
                        end
                    end
                end
                if opts.FDR,
                    if opts.designperinput,
                        for c = 1:plm.nC(1),
                            pmerged = zeros(sum(plm.Ysiz),1);
                            for y = 1:plm.nY,
                                m = y;
                                pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gtfcepperm{y}{m}{c};
                            end
                            pfdradj = fastfdr(pmerged);
                            for y = 1:plm.nY,
                                m = y;
                                palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                    sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                            end
                        end
                    else
                        for m = 1:plm.nM,
                            for c = 1:plm.nC(m),
                                pmerged = zeros(sum(plm.Ysiz),1);
                                for y = 1:plm.nY,
                                    pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gtfcepperm{y}{m}{c};
                                end
                                pfdradj = fastfdr(pmerged);
                                for y = 1:plm.nY,
                                    palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                        sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_mfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Save FWER & FDR corrected across contrasts.
    if opts.corrcon,
        fprintf('Saving p-values (corrected across contrasts).\n');
        
        % FDR correction (non-spatial stats)
        if opts.FDR,
            for y = 1:plm.nY,
                if opts.designperinput,
                    loopM = y;
                    pmerged = zeros(plm.nC(1),plm.Ysiz(y));
                else
                    loopM = 1:plm.nM;
                    pmerged = zeros(sum(plm.nC),plm.Ysiz(y));
                end
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        pmerged(j,:) = plm.Gpperm{y}{m}{c};
                        j = j + 1;
                    end
                end
                pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave(pfdradj(j,:),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_cfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        j = j + 1;
                    end
                end
            end
        end
        
        if ~ opts.accel.negbin,
            
            % FWER correction (non-spatial stats)
            for y = 1:plm.nY,
                if opts.designperinput,
                    loopM = y;
                    distmax = zeros(plm.nP{1}(1),plm.nC(1));
                else
                    loopM = 1:plm.nM;
                    distmax = zeros(plm.nP{1}(1),sum(plm.nC));
                end
                j = 1;
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gmax{y}{m}{c};
                        j = j + 1;
                    end
                end
                distmax = max(distmax,[],2);
                for m = loopM,
                    for c = 1:plm.nC(m),
                        if opts.accel.tail,
                            Ptosave = palm_pareto(plm.G{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                        elseif opts.accel.gamma,
                            Ptosave = approxgamma(plm.G{y}{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                        else
                            Ptosave = palm_datapval(plm.G{y}{m}{c},distmax,false);
                        end
                        palm_quicksave( ...
                            Ptosave,1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_cfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                end
            end
            
            % Cluster statistic
            if opts.cluster.uni.do,
                for y = 1:plm.nY,
                    if opts.designperinput,
                        loopM = y;
                        distmax = zeros(plm.nP{1}(1),plm.nC(1));
                    else
                        loopM = 1:plm.nM;
                        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
                    end
                    j = 1;
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            distmax(:,j) = plm.Gclumax{y}{m}{c};
                            j = j + 1;
                        end
                    end
                    distmax = max(distmax,[],2);
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            if opts.accel.tail,
                                Ptosave = palm_pareto(plm.Gclu{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                            elseif opts.accel.gamma,
                                Ptosave = approxgamma(plm.Gclu{y}{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                            else
                                Ptosave = palm_datapval(plm.Gclu{y}{m}{c},distmax,false);
                            end
                            palm_quicksave( ...
                                Ptosave,1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,opts.cluster.str,plm.Gname{m}{c},'_cfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                end
            end
            
            % TFCE
            if opts.tfce.uni.do,
                for y = 1:plm.nY,
                    if opts.designperinput,
                        loopM = y;
                        distmax = zeros(plm.nP{1}(1),plm.nC(1));
                    else
                        loopM = 1:plm.nM;
                        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
                    end
                    j = 1;
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            distmax(:,j) = plm.Gtfcemax{y}{m}{c};
                            j = j + 1;
                        end
                    end
                    distmax = max(distmax,[],2);
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            if opts.accel.tail,
                                Ptosave = palm_pareto(plm.Gtfce{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                            elseif opts.accel.gamma,
                                Ptosave = approxgamma(plm.Gtfce{y}{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                            else
                                Ptosave = palm_datapval(plm.Gtfce{y}{m}{c},distmax,false);
                            end
                            palm_quicksave( ...
                                Ptosave,1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_cfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                end
                if opts.FDR,
                    for y = 1:plm.nY,
                        if opts.designperinput,
                            loopM = y;
                            pmerged = zeros(plm.nC(1),plm.Ysiz(y));
                        else
                            loopM = 1:plm.nM;
                            pmerged = zeros(sum(plm.nC),plm.Ysiz(y));
                        end
                        j = 1;
                        for m = loopM,
                            for c = 1:plm.nC(m),
                                pmerged(j,:) = plm.Gtfcepperm{y}{m}{c};
                                j = j + 1;
                            end
                        end
                        pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
                        j = 1;
                        for m = loopM,
                            for c = 1:plm.nC(m),
                                palm_quicksave(pfdradj(j,:),1,opts,plm,y,m,c, ...
                                    sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_cfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                                j = j + 1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Save FWER & FDR corrected across modalities and contrasts.
    if opts.corrmod && opts.corrcon,
        fprintf('Saving p-values (corrected across modalities and contrasts).\n')
        
        % FDR correction (non-spatial stats)
        if opts.FDR,
            if opts.designperinput,
                pmerged = zeros(plm.nC(1),sum(plm.Ysiz));
            else
                pmerged = zeros(sum(plm.nC),sum(plm.Ysiz));
            end
            j = 1;
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        pmerged(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gpperm{y}{m}{c};
                        j = j + 1;
                    end
                end
            end
            pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        palm_quicksave(pfdradj(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mcfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                end
            end
        end
        
        if ~ opts.accel.negbin,
            
            % FWER correction (non-spatial stats)
            if opts.designperinput,
                distmax = zeros(plm.nP{1}(1),plm.nY*plm.nC(1));
            else
                distmax = zeros(plm.nP{1}(1),plm.nY*sum(plm.nC));
            end
            j = 1;
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        distmax(:,j) = plm.Gmax{y}{m}{c};
                        j = j + 1;
                    end
                end
            end
            distmax = max(distmax,[],2);
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        if opts.accel.tail,
                            Ptosave = palm_pareto(plm.G{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                        elseif opts.accel.gamma,
                            Ptosave = approxgamma(plm.G{y}{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                        else
                            Ptosave = palm_datapval(plm.G{y}{m}{c},distmax,false);
                        end
                        palm_quicksave( ...
                            Ptosave,1,opts,plm,y,m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},'_mcfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                    end
                end
            end
            
            % Cluster statistic
            if opts.cluster.uni.do,
                if opts.designperinput,
                    distmax = zeros(plm.nP{1}(1),plm.nY*plm.nC(1));
                else
                    distmax = zeros(plm.nP{1}(1),plm.nY*sum(plm.nC));
                end
                j = 1;
                for y = 1:plm.nY,
                    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            distmax(:,j) = plm.Gclumax{y}{m}{c};
                            j = j + 1;
                        end
                    end
                end
                distmax = max(distmax,[],2);
                for y = 1:plm.nY,
                    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            if opts.accel.tail,
                                Ptosave = palm_pareto(plm.Gclu{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                            elseif opts.accel.gamma,
                                Ptosave = approxgamma(plm.Gclu{y}{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                            else
                                Ptosave = palm_datapval(plm.Gclu{y}{m}{c},distmax,false);
                            end
                            palm_quicksave( ...
                                Ptosave,1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,opts.cluster.str,plm.Gname{m}{c},'_mcfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                end
            end
            
            % TFCE
            if opts.tfce.uni.do,
                if opts.designperinput,
                    distmax = zeros(plm.nP{1}(1),plm.nY*plm.nC(1));
                else
                    distmax = zeros(plm.nP{1}(1),plm.nY*sum(plm.nC));
                end
                j = 1;
                for y = 1:plm.nY,
                    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            distmax(:,j) = plm.Gtfcemax{y}{m}{c};
                            j = j + 1;
                        end
                    end
                end
                distmax = max(distmax,[],2);
                for y = 1:plm.nY,
                    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                    for m = loopM,
                        for c = 1:plm.nC(m),
                            if opts.accel.tail,
                                Ptosave = palm_pareto(plm.Gtfce{y}{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                            elseif opts.accel.gamma,
                                Ptosave = approxgamma(plm.Gtfce{y}{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                            else
                                Ptosave = palm_datapval(plm.Gtfce{y}{m}{c},distmax,false);
                            end
                            palm_quicksave( ...
                                Ptosave,1,opts,plm,y,m,c, ...
                                sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_mcfwep',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                        end
                    end
                end
                if opts.FDR,
                    if opts.designperinput,
                        pmerged = zeros(plm.nC(1),sum(plm.Ysiz));
                    else
                        pmerged = zeros(sum(plm.nC),sum(plm.Ysiz));
                    end
                    j = 1;
                    for y = 1:plm.nY,
                        if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                        for m = loopM,
                            for c = 1:plm.nC(m),
                                pmerged(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Gtfcepperm{y}{m}{c};
                                j = j + 1;
                            end
                        end
                    end
                    pfdradj = reshape(fastfdr(pmerged(:)),size(pmerged));
                    for y = 1:plm.nY,
                        if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                        for m = loopM,
                            for c = 1:plm.nC(m),
                                palm_quicksave(pfdradj(c,plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,y,m,c, ...
                                    sprintf('%s',opts.o,opts.tfce.str,plm.Gname{m}{c},'_mcfdrp',plm.ystr{y},plm.mstr{m},plm.cstr{m}{c}));
                            end
                        end
                    end
                end
            end
        end
    end
end

% Save NPC between modalities, corrected within contrasts
if opts.npcmod && ~ opts.npccon,
    fprintf('Saving p-values for NPC between modalities (uncorrected and corrected within contrasts).\n');
    if opts.designperinput, loopM = 1; else loopM = 1:plm.nM; end
    for m = loopM,
        if opts.designperinput, mstr = ''; else mstr = plm.mstr{m}; end
        for c = 1:plm.nC(m),
            
            % NPC p-value
            if opts.saveuncorrected,
                palm_quicksave(plm.Tpperm{m}{c},1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncp',mstr,plm.cstr{m}{c}));
            end
            
            % NPC FWER-corrected
            if opts.accel.tail,
                Ptosave = palm_pareto(plm.T{m}{c},plm.Tmax{m}{c},plm.npcrev,opts.accel.tail_thr,opts.accel.G1out);
            elseif opts.accel.gamma,
                Ptosave = approxgamma(plm.T{m}{c},plm.Tmax{m}{c},plm.npcrev,1/plm.nP{m}(c),opts.accel.G1out);
            else
                Ptosave = palm_datapval(plm.T{m}{c},plm.Tmax{m}{c},plm.npcrev);
            end
            palm_quicksave( ...
                Ptosave,1,opts,plm,[],m,c, ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fwep',mstr,plm.cstr{m}{c}));
            
            % NPC FDR
            if opts.FDR,
                palm_quicksave(fastfdr(plm.Tpperm{m}{c}),1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fdrp',mstr,plm.cstr{m}{c}));
            end
            
            % Parametric combined pvalue
            if opts.savepara && ~ plm.nonpcppara && opts.saveuncorrected,
                palm_quicksave(plm.Tppara{m}{c},1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncparap',mstr,plm.cstr{m}{c}));
            end
            
            % Cluster statistic NPC results.
            if opts.cluster.npc.do,
                
                % Cluster statistic.
                palm_quicksave(plm.Tclu{m}{c},0,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,opts.cluster.str,plm.npcstr,plm.Tname,mstr,plm.cstr{m}{c}));
                
                % Cluster statistic FWER p-value
                if opts.accel.tail,
                    Ptosave = palm_pareto(plm.Tclu{m}{c},plm.Tclumax{m}{c},false,opts.accel.tail_thr,opts.accel.G1out);
                elseif opts.accel.gamma,
                    Ptosave = approxgamma(plm.Tclu{m}{c},plm.Tclumax{m}{c},false,1/plm.nP{m}(c),opts.accel.G1out);
                else
                    Ptosave = palm_datapval(plm.Tclu{m}{c},plm.Tclumax{m}{c},false);
                end
                palm_quicksave( ...
                    Ptosave,1,opts,plm,y,m,c,...
                    sprintf('%s',opts.o,opts.cluster.str,plm.npcstr,plm.Tname,'_fwep',mstr,plm.cstr{m}{c}));
            end
            
            % TFCE NPC results.
            if opts.tfce.npc.do,
                
                % TFCE statistic.
                palm_quicksave(plm.Ttfce{m}{c},0,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,mstr,plm.cstr{m}{c}));
                
                % TFCE p-value
                if opts.saveuncorrected,
                    palm_quicksave(plm.Ttfcepperm{m}{c},1,opts,plm,[],m,c,...
                        sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_uncp',mstr,plm.cstr{m}{c}));
                end
                
                % TFCE FWER p-value
                if opts.accel.tail,
                    Ptosave = palm_pareto(plm.Ttfce{m}{c},plm.Ttfcemax{m}{c},false,opts.accel.tail_thr,opts.accel.G1out);
                elseif opts.accel.gamma,
                    Ptosave = approxgamma(plm.Ttfce{m}{c},plm.Ttfcemax{m}{c},false,1/plm.nP{m}(c),opts.accel.G1out);
                else
                    Ptosave = palm_datapval(plm.Ttfce{m}{c},plm.Ttfcemax{m}{c},false);
                end
                palm_quicksave( ...
                    Ptosave,1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_fwep',mstr,plm.cstr{m}{c}));
                
                % TFCE FDR p-value
                if opts.FDR,
                    palm_quicksave(fastfdr(plm.Ttfcepperm{m}{c}),1,opts,plm,[],m,c,...
                        sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_fdrp',mstr,plm.cstr{m}{c}));
                end
            end
        end
    end
end

% Save NPC between modalities, corrected across contrasts
if opts.npcmod && ~ opts.npccon && opts.corrcon,
    fprintf('Saving p-values for NPC between modalities (corrected across contrasts).\n');
    
    % FWER correction (non-spatial stats)
    if opts.designperinput, loopM = 1; else loopM = 1:plm.nM; end
    distmax = zeros(plm.nP{1}(1),sum(plm.nC));
    j = 1;
    for m = loopM,
        for c = 1:plm.nC(m),
            distmax(:,j) = plm.Tmax{m}{c};
            j = j + 1;
        end
    end
    distmax = max(distmax,[],2);
    for m = loopM,
        if opts.designperinput, mstr = ''; else mstr  = plm.mstr{m}; end
        for c = 1:plm.nC(m),
            if opts.accel.tail,
                Ptosave = palm_pareto(plm.T{m}{c},distmax,plm.npcrev,opts.accel.tail_thr,opts.accel.G1out);
            elseif opts.accel.gamma,
                Ptosave = approxgamma(plm.T{m}{c},distmax,plm.npcrev,1/plm.nP{1}(1),opts.accel.G1out);
            else
                Ptosave = palm_datapval(plm.T{m}{c},distmax,plm.npcrev);
            end
            palm_quicksave( ...
                Ptosave,1,opts,plm,[],m,c,...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_cfwep',mstr,plm.cstr{m}{c}));
        end
    end
    
    % FDR correction (non-spatial stats)
    if opts.FDR,
        pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                pmerged(j,:) = plm.Tpperm{m}{c};
                j = j + 1;
            end
        end
        pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave(pfdradj(j,:),1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_cfdrp',plm.mstr{m},plm.cstr{m}{c}));
                j = j + 1;
            end
        end
    end
    
    % Cluster statistic NPC
    if opts.cluster.npc.do,
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Tclumax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                if opts.accel.tail,
                    Ptosave = palm_pareto(plm.Tclu{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                elseif opts.accel.gamma,
                    Ptosave = approxgamma(plm.Tclu{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                else
                    Ptosave = palm_datapval(plm.Tclu{m}{c},distmax,false);
                end
                palm_quicksave( ...
                    Ptosave,1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,opts.cluster.str,plm.npcstr,plm.Tname,'_cfwep',plm.mstr{m},plm.cstr{m}{c}));
            end
        end
    end

    % TFCE NPC
    if opts.tfce.npc.do,
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Ttfcemax{m}{c};
                j = j + 1;
            end
        end
        distmax = max(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                if opts.accel.tail,
                    Ptosave = palm_pareto(plm.Ttfce{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                elseif opts.accel.gamma,
                    Ptosave = approxgamma(plm.Ttfce{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                else
                    Ptosave = palm_datapval(plm.Ttfce{m}{c},distmax,false);
                end
                palm_quicksave( ...
                    Ptosave,1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_cfwep',plm.mstr{m},plm.cstr{m}{c}));
            end
        end
        if opts.FDR,
            pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    pmerged(j,:) = plm.Ttfcepperm{m}{c};
                    j = j + 1;
                end
            end
            pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    palm_quicksave(pfdradj(j,:),1,opts,plm,[],m,c, ...
                        sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_cfdrp',plm.mstr{m},plm.cstr{m}{c}));
                    j = j + 1;
                end
            end
        end
    end
end

% Save NPC between contrasts, corrected within modality
if opts.npccon,
    fprintf('Saving p-values for NPC between contrasts (uncorrected and corrected within modality).\n');
    
    for j = 1:numel(plm.Tmax),
        
        % NPC p-value
        if opts.saveuncorrected,
            palm_quicksave(plm.Tpperm{j},1,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncp',plm.jstr{j}));
        end

        % NPC FWER-corrected
        if opts.accel.tail,
            Ptosave = palm_pareto(plm.T{j},plm.Tmax{j},plm.npcrev,opts.accel.tail_thr,opts.accel.G1out);
        elseif opts.accel.gamma,
            Ptosave = approxgamma(plm.T{j},plm.Tmax{j},plm.npcrev,1/plm.nP{1}(1),opts.accel.G1out);
        else
            Ptosave = palm_datapval(plm.T{j},plm.Tmax{j},plm.npcrev);
        end
        palm_quicksave( ...
            Ptosave,1,opts,plm,j,[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fwep',plm.jstr{j}));

        % NPC FDR
        if opts.FDR,
            palm_quicksave(fastfdr(plm.Tpperm{j}),1,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_fdrp',plm.jstr{j}));
        end
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara && opts.saveuncorrected,
            palm_quicksave(plm.Tppara{j},1,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,'_uncparap',plm.jstr{j}));
        end
        
        % Cluster statistic NPC results.
        if opts.cluster.npc.do,
            
            % Cluster statistic.
            palm_quicksave(plm.Tclu{j},0,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,opts.cluster.str,plm.npcstr,plm.Tname,plm.jstr{j}));
            
            % Cluster statistic FWER p-value
            if opts.accel.tail,
                Ptosave = palm_pareto(plm.Tclu{j},plm.Tclumax{j},false,opts.accel.tail_thr,opts.accel.G1out);
            elseif opts.accel.gamma,
                Ptosave = approxgamma(plm.Tclu{j},plm.Tclumax{j},false,1/plm.nP{1}(1),opts.accel.G1out);
            else
                Ptosave = palm_datapval(plm.Tclu{j},plm.Tclumax{j},false);
            end
            palm_quicksave( ...
                Ptosave,1,opts,plm,j,[],[],...
                sprintf('%s',opts.o,opts.cluster.str,plm.npcstr,plm.Tname,'_fwep',plm.jstr{j}));
        end
        
        % TFCE NPC results.
        if opts.tfce.npc.do,
            
            % TFCE statistic.
            palm_quicksave(plm.Ttfce{j},0,opts,plm,j,[],[], ...
                sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,plm.jstr{j}));
            
            % TFCE p-value
            if opts.saveuncorrected,
                palm_quicksave(plm.Ttfcepperm{j},1,opts,plm,j,[],[],...
                    sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_uncp',plm.jstr{j}));
            end
            
            % TFCE FWER p-value
            if opts.accel.tail,
                Ptosave = palm_pareto(plm.Ttfce{j},plm.Ttfcemax{j},false,opts.accel.tail_thr,opts.accel.G1out);
            elseif opts.accel.gamma,
                Ptosave = approxgamma(plm.Ttfce{j},plm.Ttfcemax{j},false,1/plm.nP{1}(1),opts.accel.G1out);
            else
                Ptosave = palm_datapval(plm.Ttfce{j},plm.Ttfcemax{j},false);
            end
            palm_quicksave( ...
                Ptosave,1,opts,plm,j,[],[],...
                sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_fwep',plm.jstr{j}));
            
            % TFCE FDR p-value
            if opts.FDR,
                palm_quicksave(fastfdr(plm.Ttfcepperm{j}),1,opts,plm,j,[],[], ...
                    sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_fdrp',plm.jstr{j}));
            end
        end
    end
end

% Save the NPC over contrasts, corrected for modalities
if ~ opts.npcmod && opts.npccon && opts.corrmod,
    fprintf('Saving p-values for NPC over contrasts (corrected across modalities).\n')
    
    % NPC FWER-corrected across modalities.
    distmax = npcextr(cat(2,plm.Tmax{:}),2);
    for y = 1:numel(plm.nY),
        if opts.accel.tail,
            Ptosave = palm_pareto(plm.T{y},distmax,plm.npcrev,opts.accel.tail_thr,opts.accel.G1out);
        elseif opts.accel.gamma,
            Ptosave = approxgamma(plm.T{y},distmax,plm.npcrev,1/plm.nP{1}(1),opts.accel.G1out);
        else
            Ptosave = palm_datapval(plm.T{y},distmax,plm.npcrev);
        end
        palm_quicksave( ...
            Ptosave,1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara && opts.saveuncorrected,
            palm_quicksave(plm.Tppara{y},1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_uncparap',plm.ystr{y}));
        end
    end
    
    % NPC FDR correction (non-spatial stats)
    if opts.FDR,
        pmerged = zeros(sum(plm.Ysiz),1);
        for y = 1:plm.nY,
            pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Tpperm{y};
        end
        pfdradj = fastfdr(pmerged);
        for y = 1:plm.nY,
            palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{y},plm.Tname,'_mfdrp',plm.ystr{y}));
        end
    end
    
    % NPC FDR-corrected across modalities.
    distmax = npcextr(cat(2,plm.Tmax{:}),2);
    for y = 1:numel(plm.nY),
        if opts.accel.tail,
            Ptosave = palm_pareto(plm.T{y},distmax,plm.npcrev,opts.accel.tail_thr,opts.accel.G1out);
        elseif opts.accel.gamma,
            Ptosave = approxgamma(plm.T{y},distmax,plm.npcrev,1/plm.nP{1}(1),opts.accel.G1out);
        else
            Ptosave = palm_datapval(plm.T{y},distmax,plm.npcrev);
        end
        palm_quicksave( ...
            Ptosave,1,opts,plm,[],[],[], ...
            sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara && opts.saveuncorrected,
            palm_quicksave(plm.Tppara{y},1,opts,plm,[],[],[], ...
                sprintf('%s',opts.o,plm.Ykindstr{y},plm.npcstr,plm.Tname,'_uncparap',plm.ystr{y}));
        end
    end
    
    % Cluster statistic NPC results.
    if opts.cluster.npc.do,
        distmax = npcextr(cat(2,plm.Tclumax{:}),2);
        for y = 1:numel(plm.nY),
            
            % Cluster statistic.
            palm_quicksave(plm.Tclu{y},0,opts,plm,y,[],[], ...
                sprintf('%s',opts.o,opts.cluster.str,plm.npcstr,plm.Tname,plm.ystr{y}));
            
            % Cluster statistic FWER p-value
            if opts.accel.tail,
                Ptosave = palm_pareto(plm.Tclu{y},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
            elseif opts.accel.gamma,
                Ptosave = approxgamma(plm.Tclu{y},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
            else
                Ptosave = palm_datapval(plm.Tclu{y},distmax,false);
            end
            palm_quicksave( ...
                Ptosave,1,opts,plm,y,[],[],...
                sprintf('%s',opts.o,opts.cluster.str,plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        end
    end
        
    % TFCE NPC results.
    if opts.tfce.npc.do,
        distmax = npcextr(cat(2,plm.Ttfce{:}),2);
        for y = 1:numel(plm.nY),
            
            % Cluster statistic.
            palm_quicksave(plm.Ttfce{y},0,opts,plm,y,[],[], ...
                sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,plm.ystr{y}));
            
            % Cluster statistic FWER p-value
            if opts.accel.tail,
                Ptosave = palm_pareto(plm.Ttfce{y},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
            elseif opts.accel.gamma,
                Ptosave = approxgamma(plm.Ttfce{y},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
            else
                Ptosave = palm_datapval(plm.Ttfce{y},distmax,false);
            end
            palm_quicksave( ...
                Ptosave,1,opts,plm,y,[],[],...
                sprintf('%s',opts.o,opts.tfce.str,plm.npcstr,plm.Tname,'_mfwep',plm.ystr{y}));
        end
        
        % NPC FDR correction TFCE
        if opts.FDR,
            pmerged = zeros(sum(plm.Ysiz),1);
            for y = 1:plm.nY,
                pmerged(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)) = plm.Ttfcepperm{y};
            end
            pfdradj = fastfdr(pmerged);
            for y = 1:plm.nY,
                palm_quicksave(pfdradj(plm.Ycumsiz(y)+1:plm.Ycumsiz(y+1)),1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,opts.tfce.str,plm.Tname,'_mfdrp',plm.ystr{y}));
            end
        end
    end
end

% Save the MV results for each contrast
if opts.MV || opts.CCA || opts.PLS,
    fprintf('Saving p-values for classical multivariate (uncorrected and corrected within contrast).\n')
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            
            % MV p-value
            if opts.saveuncorrected,
                palm_quicksave(plm.Qpperm{m}{c},1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_uncp',plm.mstr{m},plm.cstr{m}{c}));
            end
            
            % MV FDR
            if opts.FDR,
                palm_quicksave(fastfdr(plm.Qpperm{m}{c}),1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_fdrp',plm.mstr{m},plm.cstr{m}{c}));
            end
            
            % Parametric MV pvalue
            if opts.savepara && ~ plm.nomvppara && opts.saveuncorrected,
                palm_quicksave(plm.Qppara{m}{c},1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_uncparap',plm.mstr{m},plm.cstr{m}{c}));
            end
            
            % Only permutation p-value and its FDR ajustment are saved in the negative binomial mode.
            if ~ opts.accel.negbin,
                
                % MV FWER-corrected within modality and contrast.
                if opts.accel.tail,
                    Ptosave = palm_pareto(plm.Q{m}{c},plm.Qmax{m}{c},plm.mvrev{m}{c},opts.accel.tail_thr,opts.accel.G1out);
                elseif opts.accel.noperm,
                    Ptosave = [];
                elseif opts.accel.gamma,
                    Ptosave = approxgamma(plm.Q{m}{c},plm.Qmax{m}{c},plm.mvrev{m}{c},1/plm.nP{m}(c),opts.accel.G1out);
                else
                    Ptosave = palm_datapval(plm.Q{m}{c},plm.Qmax{m}{c},plm.mvrev{m}{c});
                end
                palm_quicksave( ...
                    Ptosave,1,opts,plm,[],[],[], ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_fwep',plm.mstr{m},plm.cstr{m}{c}));
                
                % Cluster statistic MV results.
                if opts.cluster.mv.do,
                    
                    % Cluster statistic.
                    palm_quicksave(plm.Qclu{m}{c},0,opts,plm,[],[],[], ...
                        sprintf('%s',opts.o,opts.cluster.str,plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{m}{c}));
                    
                    % Cluster statistic FWER p-value
                    if opts.accel.tail,
                        Ptosave = palm_pareto(plm.Qclu{m}{c},plm.Qclumax{m}{c},false,opts.accel.tail_thr,opts.accel.G1out);
                    elseif opts.accel.gamma,
                        Ptosave = approxgamma(plm.Qclu{m}{c},plm.Qclumax{m}{c},false,1/plm.nP{m}(c),opts.accel.G1out);
                    else
                        Ptosave = palm_datapval(plm.Qclu{m}{c},plm.Qclumax{m}{c},false);
                    end
                    palm_quicksave( ...
                        Ptosave,1,opts,plm,[],[],[],...
                        sprintf('%s',opts.o,opts.cluster.str,plm.mvstr,plm.Qname{m}{c},'_fwep',plm.mstr{m},plm.cstr{m}{c}));
                end
                
                % TFCE MV results.
                if opts.tfce.mv.do,
                    
                    % TFCE statistic.
                    palm_quicksave(plm.Qtfce{m}{c},0,opts,plm,[],[],[], ...
                        sprintf('%s',opts.o,opts.tfce.str,plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{m}{c}));
                    
                    % TFCE p-value
                    if opts.saveuncorrected,
                        palm_quicksave(plm.Qtfcepperm{m}{c},1,opts,plm,[],[],[],...
                            sprintf('%s',opts.o,opts.tfce.str,plm.mvstr,plm.Qname{m}{c},'_uncp',plm.mstr{m},plm.cstr{m}{c}));
                    end
                    
                    % TFCE FWER p-value
                    if opts.accel.tail,
                        Ptosave = palm_pareto(plm.Qtfce{m}{c},plm.Qtfcemax{m}{c},false,opts.accel.tail_thr,opts.accel.G1out);
                    elseif opts.accel.gamma,
                        Ptosave = approxgamma(plm.Qtfce{m}{c},plm.Qtfcemax{m}{c},false,1/plm.nP{m}(c),opts.accel.G1out);
                    else
                        Ptosave = palm_datapval(plm.Qtfce{m}{c},plm.Qtfcemax{m}{c},false);
                    end
                    palm_quicksave( ...
                        Ptosave,1,opts,plm,[],[],[], ...
                        sprintf('%s',opts.o,opts.tfce.str,plm.mvstr,plm.Qname{m}{c},'_fwep',plm.mstr{m},plm.cstr{m}{c}));
                    
                    % TFCE MV FDR
                    if opts.FDR,
                        palm_quicksave(fastfdr(plm.Qtfcepperm{m}{c}),1,opts,plm,[],[],[], ...
                            sprintf('%s',opts.o,opts.tfce.str,plm.mvstr,plm.Qname{m}{c},'_fdrp',plm.mstr{m},plm.cstr{m}{c}));
                    end
                end
            end
        end
    end
end

% Save FWER corrected across contrasts for MV.
if ( opts.MV || opts.CCA || opts.PLS) && opts.corrcon,
    fprintf('Saving p-values for classical multivariate tests (corrected across contrasts).\n')

    % FDR correction (non-spatial stats)
    if opts.FDR,
        pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                pmerged(j,:) = plm.Qpperm{m}{c};
                j = j + 1;
            end
        end
        pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                palm_quicksave(pfdradj(j,:),1,opts,plm,[],m,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_cfdrp',plm.mstr{m},plm.cstr{m}{c}));
                j = j + 1;
            end
        end
    end
    
    % Only permutation p-value and its FDR ajustment are saved in the negative binomial mode.
    if ~ opts.accel.negbin,
        
        % FWER correction (non-spatial stats)
        distmax = zeros(plm.nP{1}(1),sum(plm.nC));
        j = 1;
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                distmax(:,j) = plm.Qmax{m}{c};
                j = j + 1;
            end
        end
        if plm.mvrev{m}{c}, mvextr = @min; else mvextr = @max; end
        distmax = mvextr(distmax,[],2);
        for m = 1:plm.nM,
            for c = 1:plm.nC(m),
                if opts.accel.tail,
                    Ptosave = palm_pareto(plm.Q{m}{c},distmax,plm.mvrev{m}{c},opts.accel.tail_thr,opts.accel.G1out);
                elseif opts.accel.gamma,
                    Ptosave = approxgamma(plm.Q{m}{c},distmax,plm.mvrev{m}{c},1/plm.nP{1}(1),opts.accel.G1out);
                else
                    Ptosave = palm_datapval(plm.Q{m}{c},distmax,plm.mvrev{m}{c});
                end
                palm_quicksave( ...
                    Ptosave,1,opts,plm,[],m,c,...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},'_cfwep',plm.mstr{m},plm.cstr{m}{c}));
            end
        end
        
        % Cluster statistic MV
        if opts.cluster.mv.do,
            distmax = zeros(plm.nP{1}(1),sum(plm.nC));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    distmax(:,j) = plm.Qclumax{m}{c};
                    j = j + 1;
                end
            end
            distmax = max(distmax,[],2);
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    if opts.accel.tail,
                        Ptosave = palm_pareto(plm.Qclu{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                    elseif opts.accel.gamma,
                        Ptosave = approxgamma(plm.Qclu{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                    else
                        Ptosave = palm_datapval(plm.Qclu{m}{c},distmax,false);
                    end
                    palm_quicksave( ...
                        Ptosave,1,opts,plm,[],m,c,...
                        sprintf('%s',opts.o,opts.cluster.str,plm.mvstr,plm.Qname{m}{c},'_cfwep',plm.mstr{m},plm.cstr{m}{c}));
                end
            end
        end
        
        % TFCE MV
        if opts.tfce.mv.do,
            
            % FWER correction
            distmax = zeros(plm.nP{1}(1),sum(plm.nC));
            j = 1;
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    distmax(:,j) = plm.Qtfcemax{m}{c};
                    j = j + 1;
                end
            end
            distmax = max(distmax,[],2);
            for m = 1:plm.nM,
                for c = 1:plm.nC(m),
                    if opts.accel.tail,
                        Ptosave = palm_pareto(plm.Qtfce{m}{c},distmax,false,opts.accel.tail_thr,opts.accel.G1out);
                    elseif opts.accel.gamma,
                        Ptosave = approxgamma(plm.Qtfce{m}{c},distmax,false,1/plm.nP{1}(1),opts.accel.G1out);
                    else
                        Ptosave = palm_datapval(plm.Qtfce{m}{c},distmax,false);
                    end
                    palm_quicksave( ...
                        Ptosave,1,opts,plm,[],m,c,...
                        sprintf('%s',opts.o,opts.tfce.str,plm.mvstr,plm.Qname{m}{c},'_cfwep',plm.mstr{m},plm.cstr{m}{c}));
                end
            end
            
            % FDR correction TFCE
            if opts.FDR,
                pmerged = zeros(sum(plm.nC),plm.Ysiz(1));
                j = 1;
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        pmerged(j,:) = plm.Qtfcepperm{m}{c};
                        j = j + 1;
                    end
                end
                pfdradj = reshape(fastfdr(pmerged(:)),sum(plm.nC),plm.Ysiz(1));
                j = 1;
                for m = 1:plm.nM,
                    for c = 1:plm.nC(m),
                        palm_quicksave(pfdradj(j,:),1,opts,plm,[],m,c, ...
                            sprintf('%s',opts.o,opts.tfce.str,plm.mvstr,plm.Qname{m}{c},'_cfdrp',plm.mstr{m},plm.cstr{m}{c}));
                        j = j + 1;
                    end
                end
            end
        end
    end
end

% ==============================================================
function padj = fastfdr(pval)
% Compute FDR-adjusted p-values.
V = numel(pval);
[pval,oidx] = sort(pval);
[~,oidxR]   = sort(oidx);
padj = zeros(size(pval));
prev = 1;
for i = V:-1:1,
    padj(i) = min(prev,pval(i)*V/i);
    prev = padj(i);
end
padj = padj(oidxR);

% ==============================================================
function pvals = approxgamma(G,Gdist,rev,prepl,G1out)
% Shortcut to the calls for Gamma approximation.
if G1out,
    Gdist = Gdist(2:end,:);
end
[mu,s2,gamm1] = palm_moments(Gdist);
pvals = palm_gamma(G,mu,s2,gamm1,rev,prepl);
