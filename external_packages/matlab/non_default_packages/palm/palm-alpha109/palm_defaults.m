function opts = palm_defaults
% Set up PALM defaults.
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

% Define some defaults and organise all as a struct
opts.o                   = 'palm';             % Default output string.
opts.nP0                 = 10000;              % Number of permutations
opts.lx                  = true;               % Lexicographic permutations?
opts.cmcp                = false;              % Use Conditional Monte Carlo (ignore repeated perms)?
opts.cmcx                = false;              % Use Conditional Monte Carlo (ignore repeated elements in X)?
opts.twotail             = false;              % Do a two-tailed t-test for the t-contrasts?
opts.concordant          = false;              % Favour hypotheses with the same sign?
opts.tonly               = false;              % Run only the t-contrasts?
opts.fonly               = false;              % Run only the F-contrasts?
opts.pmethodp            = 'Guttman';          % Method to partition the model to define the permutation set
opts.pmethodr            = 'Beckmann';         % Method to partition the model for the actual regression
opts.rmethod             = 'Freedman-Lane';    % Regression/permutation method.
opts.rfallback           = 'terBraak';         % Regression/permutation method if correcting over contrasts
opts.NPC                 = false;              % Do non-parametric combination?
opts.MV                  = false;              % Do classical multivariate inference? (MANOVA/MANCOVA)
opts.CCA                 = false;              % Do canonical correlation analysis? (CCA)
opts.PLS                 = false;              % Do partial least squares regression? (PLS)
opts.mvstat              = 'auto';             % Default method for MANOVA/MANCOVA. The "auto" will use either Hotelling's T^2 or Wilks lambda depending on rank(C). For approx.noperm, this will use Pillai.
opts.npcmethod           = 'Fisher';           % Combination method.
opts.cfallback           = 'Tippett';          % (not currently used)
opts.npcmod              = false;              % FWER correction over input modalities?
opts.npccon              = false;              % FWER correction over contrasts?
opts.savepara            = false;              % Save parametric p-values too?
opts.saveglm             = false;              % Save COPEs and VARCOPEs in the 1st permutation?
opts.savecdf             = false;              % Save 1-p instead of true p.
opts.savelogp            = false;              % Convert p-values to -log10(p).
opts.savepartial         = false;              % If running NPC or MV, save partial tests?
opts.saveunivariate      = true;               % Save univariate stats when using multivariate stats?
opts.savedof             = false;              % save a file with the degrees of freedom?
opts.savemask            = false;              % Save the masks?
opts.forcemaskinter      = false;              % Force an intersection mask across modalities and evperdat?
opts.saveperms           = false;              % Save permutation images?
opts.savemetrics         = false;              % Save permutation metrics?
opts.saveuncorrected     = true;               % Save uncorrected p-values?
opts.corrmod             = false;              % FWER correction over modalities?
opts.corrcon             = false;              % FWER correction over contrasts?
opts.FDR                 = false;              % FDR adjustment?
opts.cluster.stat        = 'extent';           % Method to do the cluster-level statistic
opts.cluster.uni.do      = false;              % Do cluster statistic for the t-stat?
opts.cluster.npc.do      = false;              % Do cluster statistic for the NPC z-stat?
opts.cluster.mv.do       = false;              % Do cluster statistic for the MV z-stat?
opts.tfce.stat           = 'tfce';             % Method to do TFCE
opts.tfce.uni.do         = false;              % Do TFCE?
opts.tfce.npc.do         = false;              % Do TFCE for NPC?
opts.tfce.mv.do          = false;              % Do TFCE for MV?
opts.tfce.H              = 2;                  % TFCE H parameter
opts.tfce.E              = 0.5;                % TFCE E parameter
opts.tfce.conn           = 6;                  % TFCE connectivity neighbourhood
opts.tfce.deltah         = 'auto';             % Delta-h of the TFCE equation
opts.tableasvolume       = false;              % Has the option -tfce1d been given?
opts.useniiclass         = true;               % Use the NIFTI class (saves memory)
opts.inormal             = false;              % Do an inverse-normal transformation?
opts.inormal_meth        = 'Waerden';          % Method for the inverse-normal transformation.
opts.inormal_quanti      = true;               % Treat inputs as quantitative for the inormal?
opts.probit              = false;              % Do a probit transformation?
opts.seed                = 0;                  % Seed for the random number generator
opts.demean              = false;              % Mean-center?
opts.vgdemean            = false;              % Mean-center within VG?
opts.ev4vg               = false;              % Add one EV for each VG?
opts.removevgbysize      = 0;                  % Remove VGs smaller than a given size?
opts.zstat               = false;              % Convert G-stat to z-stat?
opts.SwE                 = false;              % Whole block shuffling (for legacy format)?
opts.pearson             = false;              % Compute the Pearson's correlation coefficient (R^2 if rank(C)>1)?
opts.noranktest          = false;              % Don't test the rank(Y) before doing MANOVA/MANCOVA.
opts.evperdat            = false;              % Use one (single!) EV per datum?
opts.transposedata       = false;              % transpose data if 2D?
opts.verbosefilenames    = false;              % use filenames with _i%d, _c%d, etc when there are more than one?
opts.syncperms           = false;              % synchronize permutations (this only affects designs that can actually be synced
opts.designperinput      = false;              % use one design for each input?
opts.showprogress        = true;               % print progress as the permutations are performed? (use -quiet to disable)
opts.inputmv             = false;              % treat each column of the sole input as a separate input for MV/NPC/CCA?
opts.reversemasks        = false;              % Reverse the masks?
opts.precision           = [];                 % Precision? Can be 'single', 'double', or [] for what the file defines.
% Approximation options:
opts.accel.negbin        = 0;                  % Run a negbin scheme
opts.accel.negbin_nexced = 2;                  % Number of exceedances for the negbin mode
opts.accel.tail          = false;              % Use a Pareto tail approximation to the FWER-adjusted p-values?
opts.accel.tail_thr      = 0.10;               % P-values below this will be approximated using tail approximation
opts.accel.noperm        = false;              % Compute approximate p-vals without doing any permutation
opts.accel.gamma         = false;              % Do a gamma-fit after just a few permutations
opts.accel.lowrank       = false;              % Try a low rank approximation
opts.accel.lowrank_val   = NaN;                % Use NaN for N*(N+1)/2. Values <=1 are fraction of voxels to use; vals >1 are number of voxels to use.
opts.accel.lowrank_recon = false;              % Reconstruct past permutations in the new basis? This is very slow.
opts.accel.G1out         = false;              % Exclude (true) or not (false) the unpermuted statistic in the null distribution for tail and gamma?

% Missing data options:
opts.missingdata         = false;              % Are there missing data?
opts.mcar                = false;              % Data missing completely at random?
opts.npcmethodmiss       = 'Fisher';           % Combination method for missing data.

% Note that there are no adjustable defaults for EE, ISE, whole or within.
% These are hard coded and not meant to be ever changed (EE is default, within-block is also default).
