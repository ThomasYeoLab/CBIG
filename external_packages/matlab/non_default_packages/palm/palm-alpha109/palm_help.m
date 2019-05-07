function varargout = palm_help(varargin)
% Shows a help text. Call palm_help('logo') to show just the logo.
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Mar/2014
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

if nargin == 0 || strcmpi(varargin{1},'basic'),
    showlogo;
    basic_help;
    showdate;
elseif strcmpi(varargin{1},'advanced'),
    showlogo;
    advanced_help;
    showdate;
elseif strcmpi(varargin{1},'logo'),
    showlogo;
elseif strcmpi(varargin{1},'date'),
    showdate;
elseif strcmpi(varargin{1},'version'),
    if nargout == 0,
        fprintf(showversion);
    else
        varargout{1} = showversion;
    end
end

% ==============================================================
function basic_help
% Show the most common options.

fprintf('\nThe main options are:\n\n');

fprintf('-i <file> : Input(s). More than one can be specified, each one preceded\n');
fprintf('	by its own -i. All input files must contain the same number of\n');
fprintf('	observations (e.g., the same number of subjects). Except for NPC\n');
fprintf('	and MV, mixing is allowed (e.g., voxelwise, vertexwise and\n');
fprintf('	non-imaging data can be all loaded at once, and later will be all\n');
fprintf('	corrected across) if the option "-corrmod" is used.\n\n');

fprintf('-m <file> : Mask(s). Either one for all inputs, or one per input,\n');
fprintf('	supplied in the same order as the respective -i appear.\n\n');

fprintf('-s <filesurf> [filearea] : The first argument is the surface file\n');
fprintf('	itself. The second is an optional area-per-vertex or area-per-\n');
fprintf('	face file, or simply a number. If only the surface file is\n');
fprintf('	provided, its area is calculated and used for the computation\n');
fprintf('	of spatial statistics (cluster extent and TFCE). If the second\n');
fprintf('	argument is given, it should contain the areas, which are then\n');
fprintf('	used (e.g., average areas from native geometry after areal\n');
fprintf('	interpolation). Alternatively, if the areas are not meaningful\n');
fprintf('	for cluster extent or TFCE, this argument can be simply a number,\n');
fprintf('	such as "1", which is then used as the area of all vertices or\n');
fprintf('	faces.\n\n');

fprintf('-d <file> : Design matrix. It can be in csv format, or in FSL''s vest\n');
fprintf('	format. For information on how to construct the design matrix,\n');
fprintf('	see the FSL GLM manual.\n\n');

fprintf('-t <file> : t-contrasts file, in csv or vest format (the format used by\n');
fprintf('	FSL). The option -t can be used more than once, so that more than\n');
fprintf('	one t-contrasts file can be loaded.\n\n');

fprintf('-f <file> : F-contrasts file, in csv or vest format. The option -f can\n');
fprintf('	be used more than once, so that more than one F-contrasts file\n');
fprintf('	can be loaded. Each file supplied with a -f corresponds to the\n');
fprintf('	file supplied with the option -t in the same order. The option\n');
fprintf('	-f cannot be used more than the number of times the option\n');
fprintf('	-t was used.\n\n');

fprintf('-fonly : Run only the F-contrasts, not the t-contrasts.\n\n');

fprintf('-n <integer> : Number of permutations. Use -n 0 to run all permutations\n');
fprintf('	and/or sign-flips exhaustively. Default is 10000.\n\n');

fprintf('-eb <file> : Exchangeability blocks file, in csv or vest format. If\n');
fprintf('	omitted, all observations are treated as exchangeable and\n');
fprintf('	belonging to a single large exchangeability block.\n\n');

fprintf('-within : If the file supplied with -eb has a single column, this option\n');
fprintf('	runs within-block permutation (default). Can be used with "-whole".\n\n');

fprintf('-whole : If the file supplied with -eb has a single column, this option\n');
fprintf('	runs whole-block permutation. Can be used with "-within".\n\n');

fprintf('-ee : Assume exchangeable errors (EE), to allow permutations.\n\n');

fprintf('-ise : Assume independent and symmetric errors (ISE), to allow\n');
fprintf('	sign-flipping.\n\n');

fprintf('-vg <file> : Variance groups file, in csv or vest format. If omitted,\n');
fprintf('	all observations are assumed to belong to the same variance group\n');
fprintf('	(i.e., the data is treated as homoscedastic. Use "-vg auto" to\n');
fprintf('	define automatically using a structure that is compatible with the\n');
fprintf('	exchangeability blocks (option -eb).\n\n');

fprintf('-npcmethod <method> : Specify the combining function for NPC (Non-\n');
fprintf('	Parametric Combination), which can be one of: Tippett, Fisher,\n');
fprintf('	Pearson-David, Stouffer, Wilkinson <alpha>, Winer, Edgington,\n');
fprintf('	Mudholkar-George, Friston, Darlington-Hayes <r>, Zaykin <alpha>,\n');
fprintf('	Dudbridge-Koeleman <r>, Dudbridge-Koeleman2 <r> <alpha>,\n');
fprintf('	Taylor-Tibshirani or Jiang <alpha>. Default is Fisher. Note that\n');
fprintf('	some methods require 1 or 2 additional parameters to be provided.\n');
fprintf('	All methods except Darlington-Hayes and Jiang can also be used to\n');
fprintf('	produce parametric p-values and spatial statistics.\n\n');

fprintf('-npcmod : Enable NPC over modalities.\n\n');

fprintf('-npccon : Enable NPC over contrasts.\n\n');

fprintf('-npc : Shortcut to "-npcmethod <default method> -npcmod".\n\n');

fprintf('-mv <statistic> <k>: Do classical multivariate analysis (MV), such as\n');
fprintf('	MANOVA and MANCOVA, using the the specified statistic, which can\n');
fprintf('	be one of: Wilks, HotellingTsq, Lawley, Pillai, Roy_ii, Roy_iii,\n');
fprintf('	CCA, or PLS. All but Roy_iii can be used with spatial statistics.\n');
fprintf('	Alternatively, use CCA to do a Canonical Correlation Analysis, or\n');
fprintf('	PLS to do Partial Least Squares regression, with a permutation test\n');
fprintf('	on the indicated k-th canonical correlation (default is 1).\n\n');

fprintf('-pearson : Instead of t, F, v or G, compute either the Pearson"s\n');
fprintf('	correlation coefficient, r (if the constrast has rank = 1), or the\n');
fprintf('	coefficient of determination R^2 (if the constrast has rank > 1).\n');
fprintf('	For the contrasts in which some EVs are zeroed out, this option\n');
fprintf('	computes the multiple correlation coefficient (or R^2)\n');
fprintf('	corresponding to the EVs of interest.\n\n');

fprintf('-T : Enable TFCE inference for univariate (partial) tests, as well as\n');
fprintf('	for NPC and/or MV if these options have been enabled.\n\n');

fprintf('-C <real> : Enable cluster inference for univariate (partial) tests,\n');
fprintf('	with the supplied cluster-forming threshold (supplied as the\n');
fprintf('	equivalent z-score), as well as for NPC and/or MV if these options\n');
fprintf('	have been enabled. Use preferably values >3.\n\n');

fprintf('-Cstat <name> : Choose which cluster statistic should be used. Accepted\n');
fprintf('	statistics are "extent" and "mass" (see the source code for\n');
fprintf('	experimental possibilities).\n\n');

fprintf('-tfce1D : Set TFCE parameters for 1D data (synchronised timeseries)\n');
fprintf('	i.e., H = 2, E = 2, C = 6. Use this option together with -T.\n\n');

fprintf('-tfce2D : Set TFCE parameters for 2D data (surface, TBSS)  i.e.,\n');
fprintf('	H = 2, E = 1, C = 26. Use this option together with -T.\n\n');

fprintf('-corrmod : Apply FWER-correction of p-values over multiple modalities.\n\n');

fprintf('-corrcon : Apply FWER-correction of p-values over multiple contrasts.\n');
fprintf('	Because multivariate and F-stats can have different df, this option\n');
fprintf('	automatically enables "-zstat".\n\n');

fprintf('-fdr : Produce FDR-adjusted p-values.\n\n');

fprintf('-o <string> : Output prefix. It may itself be prefixed by a path.\n');
fprintf('	Default is palm.\n\n');

fprintf('-save1-p : Save (1-p) instead of the actual p-values.\n\n');

fprintf('-logp : Save the output p-values as -log(p) (or -log(1-p) if the\n');
fprintf('	option -save1-p is used).\n\n');

fprintf('-demean : Mean center the data, as well as all columns of the design\n');
fprintf('	matrix. If the design has an intercept, the intercept is removed.\n\n');

fprintf('-twotail : Run two-tailed tests for all the t-contrasts instead of\n');
fprintf('	one-tailed. If NPC is used, it naturally becomes two-tailed.\n\n');

fprintf('-concordant : For the NPC, favour alternative hypotheses with\n');
fprintf('	concordant signs. Cannot be used with "-twotail".\n\n');

fprintf('-approx/-accel <method> : Run one of various acceleration methods\n');
fprintf('	available. Legal options for <method> are below. Some may:\n');
fprintf('	accept extra parameters in brackets [ ], that are optional:\n');
fprintf('	"negbin [n]"        : Negative binomial, with [n] exceedances\n');
fprintf('	                      (default n = 2).\n');
fprintf('	"tail [pthr] [out]" : Tail approximation using a GPD for p-values\n');
fprintf('	                      below pthr (default pthr = 0.10), and excluding\n');
fprintf('	                      or not the unpermuted case in the null\n');
fprintf('	                      distribution (default out = false).\n');
fprintf('	"noperm"            : Compute permutation p-values without\n');
fprintf('	                      permutations.\n');
fprintf('	"gamma [out]"       : Gamma distribution approximation, excluding\n');
fprintf('	                      or not the unpermuted case in the null\n');
fprintf('	                      distribution (default out = false).\n');
fprintf('	"lowrank [val]"     : Low rank matrix completion: subsample a\n');
fprintf('	                      fraction(if val <= 1) or an absolute number\n');
fprintf('	                      of tests (val > 1) that undergo full testing.\n');
fprintf('	                      The remainder are computed via low rank matrix\n');
fprintf('	                      filling. Default val = NaN, such that it''s\n');
fprintf('	                      computed as N*(N+1)/2, where N is the\n');
fprintf('	                      number of subjects.\n');
fprintf('	If unsure about which option to use, and if unwilling to read the\n');
fprintf('	full paper, go with "-accel tail", setting a relatively small number\n');
fprintf('	of permutations, such as "-n 500". If using "-nouncorrected", or\n');
fprintf('	if the images are small (or for non-images), this number can\n');
fprintf('	be increased without hitting memory limits.\n\n');

fprintf('-reversemasks : Reverse 0/1 in the masks, so that the zero values\n');
fprintf('	are then used to select the voxels/vertices/faces.\n\n');

fprintf('-quiet : Don''t show progress as the shufflings are performed.\n\n');

fprintf('-advanced : Show advanced options.\n\n');


% ==============================================================
function advanced_help
% Show advanced options.

fprintf('\nThe advanced or less commonly used options are:\n\n');

fprintf('-imiss <file> : Missing data indicators for the input(s).\n\n');

fprintf('-dmiss <file> : Missing data indicators for the design(s).\n\n');

fprintf('-con <file1> <file2> : Contrast file(s) in .mset format. For hypotheses\n');
fprintf('	of the form H0: C''*Psi*D, file1 contains a set of C contrasts, and\n');
fprintf('	file2 (optional) contains a set of D contrasts.\n\n')

fprintf('-tonly : Run only the t-contrasts, not the F-contrasts.\n\n');

fprintf('-cmcp : Ignore the possibility of repeated permutations. This option\n');
fprintf('	will be ignored if the number of shufflings chosen is larger than the\n');
fprintf('	maximum number of possible shufflings, in which case all possible\n');
fprintf('	shufflings will be performed.\n\n');

fprintf('-cmcx : Ignore the possibility of repeated rows in X when\n');
fprintf('	constructing the set of permutations, such that each row is\n');
fprintf('	treated as unique, regardless of potential repetitions (ties).\n\n');

fprintf('-conskipcount <integer> : Normally the contrasts are numbered from 1, but\n');
fprintf('	this option allows staring the counter from the specified number.\n');
fprintf('	This option doesn"t affect which contrasts are performed.\n\n');

fprintf('-Cuni <real> : Enable cluster statistics for t contrasts for univariate\n');
fprintf('	(partial) tests, with the supplied cluster-forming threshold (as\n');
fprintf('	a z-score).\n\n');

fprintf('-Cnpc <real> : Enable cluster statistics for t contrasts for NPC, with the\n');
fprintf('	supplied cluster-forming threshold (as a z-score).\n\n');

fprintf('-Cmv <real> : Enable cluster statistics for t contrasts for MV, with the\n');
fprintf('	supplied cluster-forming threshold (as a z-score).\n\n');

fprintf('-designperinput : Use one design file for each input modality.\n\n');

fprintf('-ev4vg : Add to the design one EV modelling the mean of each VG.\n\n');

fprintf('-evperdat : Treat the design matrix supplied with -d as having one column\n');
fprintf('	for each column in the observed data (entered with -i). This\n');
fprintf('	enables voxelwise/facewise/vertexwise regressors.\n\n');

fprintf('-forceintersectionmask : Force the use of an intersection mask across\n');
fprintf('	inputs and pointwise EVs, even if not using NPC or MV.\n\n');

fprintf('-inormal : Apply an inverse-normal transformation to the data.\n');
fprintf('	This procedure can go faster if the data is known to be quantitative\n');
fprintf('	(in which case, use "-inormal quanti"). There are four different\n');
fprintf('	methods available, which can be specified as "-inormal <method>" or\n');
fprintf('	"-inormal quanti <method>". The methods are "Waerden" (default),\n');
fprintf('	"Blom", "Tukey" and "Bliss".\n\n');

fprintf('-probit : Apply a probit transformation to the data.\n\n');

fprintf('-inputmv : Treat the (sole) input as multivariate, that is, each column\n');
fprintf('	is a variable in a multivariate model, as opposed to independent\n');
fprintf('	univariate tests. Useful with non-imaging data. When used, the\n');
fprintf('	option "-nounivariate" is automatically set.\n\n');

fprintf('-mcar : Assume that data is missing completely at random.\n\n');

fprintf('-noniiclass : Do not use the NIFTI class (use this option only if you\n');
fprintf('	have problems and even so, only for small datasets).\n\n');

fprintf('-precision <string> : Precision ("single" or "double") for input files.\n\n');

fprintf('-noranktest : For MV, don''t check the rank of the data before trying to\n');
fprintf('	compute the multivariate statistics.\n\n');

fprintf('-nounivariate : Don''t save univariate resuts.\n\n');

fprintf('-nouncorrected : Don''t save uncorrected p-values.\n\n');

fprintf('-saveuncorrected : Save uncorrected p-values if an earlier option disabled it.\n\n');

fprintf('-pmethodp : Partition method used when defining the set of permutations.\n');
fprintf('	Cab be "Guttman", "Beckmann", "Ridgway" or "None".\n');
fprintf('	Default is "Beckmann"\n\n');

fprintf('-pmethodr : Partition method used during the regression. Valid values\n');
fprintf('	are "Guttman", "Beckmann", "Ridgway" or "None".\n');
fprintf('	Default is "Beckmann"\n\n');

fprintf('-removevgbysize <integer> : Remove from the data and design those\n');
fprintf('	observations that are in VGs of size smaller or equal than specified.\n');
fprintf('	This is specially useful with the option "-vg auto".\n\n')

fprintf('-rmethod <string> : Method for regression/permutation. It can be one of\n');
fprintf('	Freedman-Lane, Dekker, terBraak, Manly, Draper-Stoneman,\n');
fprintf('	Still-White and Huh-Jhun. Default, and recommended, is\n');
fprintf('	Freedman-Lane.\n\n');

fprintf('-savedof : Save file with degrees of freedom.\n\n');

fprintf('-savemask : Save the effective masks used for each modality, as well as\n');
fprintf('	an intersection mask used for NPC and/or MV of these have been\n');
fprintf('	selected.\n\n');

fprintf('-savemetrics : Save a csv file with various permutation metrics.\n\n');

fprintf('-saveparametric : Save also uncorrected parametric p-values. These are only\n');
fprintf('	valid if all assumptions are met, including iid and normality.\n\n');

fprintf('-saveglm : Save COPEs, VARCOPEs, and Cohen''s d in the first permutation\n');
fprintf('	for the contrasts that have rank = 1.\n\n');

fprintf('-saveperms : Save one statistic image per permutation, as well as a csv\n');
fprintf('	file with the permutation indices (one permutation per row, one\n');
fprintf('	index per column; sign-flips are represented by the negative\n');
fprintf('	indices). Use cautiously, as this may consume large amounts of\n');
fprintf('	disk space.\n\n');

fprintf('-seed <positive> : Seed used for the random number generator. Use a\n');
fprintf('	positive integer up to 2^32. Default is 0. To start with a random\n');
fprintf('	number, use the word "twist" instead of the integer. Note that a\n');
fprintf('	given seed used in Matlab isn''t equivalent to the same seed used\n');
fprintf('	in Octave.\n\n');

fprintf('-syncperms : If possible, use synchronized permutations even they wouldn''t\n');
fprintf('	be used by default.\n\n');

fprintf('-subjidx <file> : Indices of input data to keep in the design.\n\n');

fprintf('-Tuni : Enable TFCE inference for univariate (partial) tests.\n\n');

fprintf('-Tnpc : Enable TFCE inference for NPC.\n\n');

fprintf('-Tmv : Enable TFCE inference for MV.\n\n');

fprintf('-tfce_H <real> : Set the TFCE H parameter (height power).\n\n');

fprintf('-tfce_E <real> : Set the TFCE E parameter (extent power).\n\n');

fprintf('-tfce_C <integer> : Set the TFCE C parameter (connectivity).\n\n');

fprintf('-tfce_dh <real> : Set the "delta h" for the calculation of TFCE.\n');
fprintf('	Default is "auto".\n\n');

fprintf('-tableasvolume : Treat tables (e.g., CSV inputs) as volumes, such\n');
fprintf('	that TFCE can be calculated. This is useful for TFCE over timeseries.\n\n');

fprintf('-transposedata : For input data (-i) that are csv tables (2D), transpose\n');
fprintf('	rows and columns.\n\n');

fprintf('-verbosefilenames : Use lengthy filenames, even if the numbering go up\n');
fprintf('	to 1 only.\n\n');

fprintf('-vgdemean : Mean center the data, as well as all columns of the design\n');
fprintf('	matrix, within each VG. Intercepts are removed.\n\n');

fprintf('-zstat : Convert the original statistic (t, F, v, G, r, R2, or any of\n');
fprintf('	the MV statistics) to a normally distributed, z-statistic.\n\n');


% ==============================================================
function showlogo
% Show the logo
fprintf('=======================================================================\n');
fprintf('             ___         ___                         ___\n');
fprintf('            /  /\\       /  /\\                       /__/\\\n');
fprintf('           /  /::\\     /  /::\\                     |  |::\\\n');
fprintf('          /  /:/\\:\\   /  /:/\\:\\    ___     ___     |  |:|:\\\n');
fprintf('         /  /:/~/:/  /  /:/~/::\\  /__/\\   /  /\\  __|__|:|\\:\\\n');
fprintf('        /__/:/ /:/  /__/:/ /:/\\:\\ \\  \\:\\ /  /:/ /__/::::| \\:\\\n');
fprintf('        \\  \\:\\/:/   \\  \\:\\/:/__\\/  \\  \\:\\  /:/  \\  \\:\\~~\\__\\/\n');
fprintf('         \\  \\::/     \\  \\::/        \\  \\:\\/:/    \\  \\:\\\n');
fprintf('          \\  \\:\\      \\  \\:\\         \\  \\::/      \\  \\:\\\n');
fprintf('           \\  \\:\\      \\  \\:\\         \\__\\/        \\  \\:\\\n');
fprintf('            \\__\\/       \\__\\/                       \\__\\/\n');
fprintf('\n');
fprintf('=======================================================================\n');
fprintf('                 Permutation Analysis of Linear Models\n');
fprintf('=======================================================================\n');


% ==============================================================
function showdate
% Show the date
fprintf('_____________________________________\n');
fprintf('Anderson M. Winkler\n');
fprintf('FMRIB / University of Oxford\n');
fprintf('%s',showversion);
fprintf('http://www.fmrib.ox.ac.uk/fsl\n');

function vstr = showversion
% Read the file with the version
fid = fopen(fullfile(fileparts(mfilename('fullpath')),'palm_version.txt'),'r');
vstr = textscan(fid,'%s');
fclose(fid);

% Assemble back as a string
vstr = sprintf('%s ',vstr{1}{:});
vstr = sprintf('%s\n',vstr(1:end-1));
