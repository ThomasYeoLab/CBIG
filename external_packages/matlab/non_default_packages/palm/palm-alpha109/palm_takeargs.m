function [opts,plm] = palm_takeargs(varargin)
% Handle the inputs for PALM.
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

% Load the defaults
opts = palm_defaults;

% As varargin is actually from another function, fix it.
if nargin == 1,
    if exist(varargin{1},'file'),
        vararginx = palm_configrw(varargin{1});
    else
        error('Unknown option or file not found: %s',varargin{1});
    end
else
    vararginx = varargin;
    idxa = find(strcmpi(vararginx,'-o'));
    if isempty(idxa),
        otmp = opts.o;
    else
        otmp = vararginx{idxa+1};
    end
    if ~ strcmp(otmp(end),'_'),
        otmp = horzcat(otmp,'_');
    end
    cfgname = horzcat(otmp,'palmconfig.txt');
    [opth,~,~] = fileparts(cfgname);
    if ~isempty(opth) && ~exist(opth,'dir'),
        mkdir(opth);
    end
    palm_configrw(vararginx,cfgname);
end

% Number of input images/masks/surfaces
% These are NOT meant to be edited.
Ni     = sum(strcmp(vararginx,'-i'));        % number of data inputs
Nm     = sum(strcmp(vararginx,'-m'));        % number of masks
Ns     = sum(strcmp(vararginx,'-s'));        % number of surfaces
Nd     = sum(strcmp(vararginx,'-d'));        % number of design files
Nimiss = sum(strcmp(vararginx,'-imiss'));    % number of missing indicators for inputs
Ndmiss = sum(strcmp(vararginx,'-dmiss'));    % number of missing indicators for designs
Nt     = sum(strcmp(vararginx,'-t'));        % number of t-contrast files
Nf     = sum(strcmp(vararginx,'-f'));        % number of F-test files
Ncon   = sum(strcmp(vararginx,'-con'));      % number of contrast files (t or F, mset format)
Nevd   = sum(strcmp(vararginx,'-evperdat')); % number of EV per datum inputs
opts.i     = cell(Ni,1);   % Input files (to constitute Y later)
opts.m     = cell(Nm,1);   % Mask file(s)
opts.s     = cell(Ns,1);   % Surface file(s)
opts.sa    = cell(Ns,1);   % Area file(s) or weight(s)
opts.d     = cell(Nd,1);   % Design file(s)
opts.imiss = cell(Nd,1);   % Design file(s)
opts.dmiss = cell(Nd,1);   % Design file(s)
opts.t     = cell(Nt,1);   % t contrast file(s)
opts.f     = opts.t;       % F contrast file(s)
opts.Ccon  = cell(Ncon,1); % Contrast file(s) (t or F, mset format)
opts.Dcon  = cell(Ncon,1); % Contrast file(s) (multivariate, mset format)
opts.eb       = [];       % File with definition of exchangeability blocks
opts.vg       = [];       % File with definition of variance groups
opts.EE       = false;    % To be filled below (don't edit this!)
opts.ISE      = false;    % To be filled below (don't edit this!)
opts.within   = false;    % To be filled below (don't edit this!)
opts.whole    = false;    % To be filled below (don't edit this!)
opts.conskipcount = 0;    % When saving the contrasts, skip how many from 1?
opts.singlevg = true;     % Make sure that sigle VG will be used if nothing is supplied (this is NOT a "default" setting, and it's not a setting at all, but hard coded. Don't edit it!)
opts.subjidx  = [];       % Filename of the indices of subjects to keep
plm.subjidx   = [];       % Indices of subjects to keep

% These are to be incremented below
i = 1; m = 1; d = 1;
t = 1; f = 1; s = 1;
con = 1; ev = 1;
imiss = 1; dmiss = 1;

% Remove trailing empty arguments. This is useful for some Octave versions.
while numel(vararginx) > 0 && isempty(vararginx{1}),
    vararginx(1) = [];
end
narginx = numel(vararginx);

% Take the input arguments
a = 1;
while a <= narginx,
    switch vararginx{a},
        case {'-help','-?','-basic','-advanced'},
            
            % Do nothing, as these options are parsed separately,
            % and should anyway be given without any other argument.
            a = a + 1;
            
        case '-i', % basic
            
            % Get the filenames for the data.
            opts.i{i} = vararginx{a+1};
            i = i + 1;
            a = a + 2;
            
        case '-m', % basic
            
            % Get the filenames for the masks, if any.
            opts.m{m} = vararginx{a+1};
            m = m + 1;
            a = a + 2;

        case {'-s','-surf'}, % basic
            
            % Get the filenames for the surfaces, if any.
            opts.s{s}  = vararginx{a+1};
            if nargin == a+1 || (narginx>a+1 && strcmp(vararginx{a+2}(1),'-')),
                opts.sa{s} = [];
                a = a + 2;
            else
                opts.sa{s} = vararginx{a+2};
                a = a + 3;
            end
            s = s + 1;
            
        case '-d', % basic
            
            % Get the design matrix file.
            opts.d{d} = vararginx{a+1};
            d = d + 1;
            a = a + 2;
            
        case '-evperdat', % advanced
            
            % Use one EV per datum?
            opts.evperdat = true;
            opts.evdatfile{ev} = vararginx{a+1};
            if nargin == a + 1 || ...
                    ischar(vararginx{a+2}) && ...
                    strcmpi(vararginx{a+2}(1),'-'),
                opts.evpos{ev}(1) = 1;  % EV position
                opts.evpos{ev}(2) = 1;  % Design number
                a = a + 2;
            elseif nargin == a + 2 || ...
                    ischar(vararginx{a+3}) && ...
                    strcmpi(vararginx{a+3}(1),'-'),
                if ischar(vararginx{a+2}), % EV position
                    opts.evpos{ev}(1) = eval(vararginx{a+2});
                else
                    opts.evpos{ev}(1) = vararginx{a+2};
                end
                opts.evpos{ev}(2) = 1;  % Design number
                a = a + 3;
            else
                if ischar(vararginx{a+2}), % EV position
                    opts.evpos{ev}(1) = eval(vararginx{a+2});
                else
                    opts.evpos{ev}(1) = vararginx{a+2};
                end
                if ischar(vararginx{a+3}), % Design number
                    opts.evpos{ev}(2) = eval(vararginx{a+3});
                else
                    opts.evpos{ev}(2) = vararginx{a+3};
                end
                a = a + 4;
            end
            ev = ev + 1;
            
        case '-imiss', % basic
            
            % Get the filenames for the missing data indicators (inputs).
            opts.imiss{imiss} = vararginx{a+1};
            imiss = imiss + 1;
            a = a + 2;
            
        case '-dmiss', % basic
            
            % Get the filenames for the missing data indicators (designs).
            opts.dmiss{dmiss} = vararginx{a+1};
            dmiss = dmiss + 1;
            a = a + 2;
            
        case '-mcar',
            
            % For the missing data, treat as missing completely at random.
            opts.mcar = true;
            a = a + 1;
            
        case '-t', % basic
            
            % Get the t contrast files.
            opts.t{t} = vararginx{a+1};
            t = t + 1;
            a = a + 2;
            
        case '-f', % basic
            
            % Get the F contrast files.
            if t == 1,
                error('The option "-f" cannot be specified before its respective "-t".');
            end
            opts.f{t-1} = vararginx{a+1};
            a = a + 2;
            
        case '-con', % advanced
            
            % Get the contrast files from an .mset file or
            % pair of files. If a pair, the 1st is for Cset
            % and the second for Dset.
            opts.Ccon{con} = vararginx{a+1};
            if nargin == a + 1 || ...
                    ischar(vararginx{a+2}) && ...
                    strcmpi(vararginx{a+2}(1),'-'),
                opts.Dcon{con} = [];
                a = a + 2;
            else
                opts.Dcon{con} = vararginx{a+2};
                a = a + 3;
            end
            con = con + 1;
            
        case '-conskipcount', % advanced
            
            % Numbers to skip when saving the contrasts
            opts.conskipcount = vararginx{a+1};
            if ischar(opts.conskipcount),
                opts.conskipcount = str2double(opts.conskipcount);
            end
            a = a + 2;
            
        case '-tonly', % advanced
            
            % Run only the t-contrasts
            opts.tonly = true;
            a = a + 1;
            
        case '-fonly', % basic
            
            % Run only the F-contrasts
            opts.fonly = true;
            a = a + 1;
            
        case '-eb', % basic
            
            % Get the exchangeability blocks file.
            opts.eb = vararginx{a+1};
            a = a + 2;
            
        case '-vg' % basic
            
            % Get the variance groups file.
            opts.vg = vararginx{a+1};
            if     ischar(opts.vg) && ...
                    any(strcmpi(opts.vg,{'single'})),
                opts.vg = 'single';
                opts.singlevg = true;
            elseif ischar(opts.vg) && ...
                    any(strcmpi(opts.vg,{'auto','automatic'})),
                opts.vg = 'auto';
                opts.singlevg = false;
            else
                opts.singlevg = false;
            end
            a = a + 2;
            
        case '-swe', % advanced
            
            % Compute one (of various possible) sandwich estimators
            opts.SwE = true;
            a = a + 1;
            
        case '-o', % basic
            
            % Output prefix for the files to be saved.
            opts.o = vararginx{a+1};
            a = a + 2;
            
        case '-n', % basic
            
            % Number of permutations
            opts.nP0 = vararginx{a+1};
            if ischar(opts.nP0),
                opts.nP0 = str2double(opts.nP0);
            end
            a = a + 2;
            
        case '-C', % basic
            
            % Threshold for cluster extent, univariate, NPC and MV
            opts.cluster.uni.do = true;
            opts.cluster.uni.thr = vararginx{a+1};
            if ischar(opts.cluster.uni.thr),
                opts.cluster.uni.thr = str2double(opts.cluster.uni.thr);
            end
            opts.cluster.npc.do = true;
            opts.cluster.npc.thr = vararginx{a+1};
            if ischar(opts.cluster.npc.thr),
                opts.cluster.npc.thr = str2double(opts.cluster.npc.thr);
            end
            opts.cluster.mv.do  = true;
            opts.cluster.mv.thr = vararginx{a+1};
            if ischar(opts.cluster.mv.thr),
                opts.cluster.mv.thr = str2double(opts.cluster.mv.thr);
            end
            a = a + 2;
            
        case '-Cuni', % advanced
            
            % Threshold for cluster statistic, univariate
            opts.cluster.uni.do = true;
            opts.cluster.uni.thr = vararginx{a+1};
            if ischar(opts.cluster.uni.thr),
                opts.cluster.uni.thr = str2double(opts.cluster.uni.thr);
            end
            a = a + 2;
            
        case '-Cnpc', % advanced
            
            % Threshold for cluster statistic, NPC
            opts.NPC = true;
            opts.cluster.npc.do = true;
            opts.cluster.npc.thr = vararginx{a+1};
            if ischar(opts.cluster.npc.thr),
                opts.cluster.npc.thr = str2double(opts.cluster.npc.thr);
            end
            a = a + 2;
            
        case '-Cmv', % advanced
            
            % Threshold for cluster statistic, MV
            opts.MV = true;
            opts.cluster.mv.do = true;
            opts.cluster.mv.thr = vararginx{a+1};
            if ischar(opts.cluster.mv.thr),
                opts.cluster.mv.thr = str2double(opts.cluster.mv.thr);
            end
            a = a + 2;
            
        case '-Cstat', % advanced
            
            % Type of cluster statistic
            opts.cluster.stat = vararginx{a+1};
            if ~ any(strcmp(opts.cluster.stat,{'extent','mass','density','tippett','pivotal'})),
                error('Cluster statistic "%s" unknown.',opts.cluster.stat);
            end
            a = a + 2;
            
        case '-T', % basic
            
            % Do TFCE for univariate, NPC and MV?
            opts.tfce.uni.do = true;
            opts.tfce.npc.do = true;
            opts.tfce.mv.do = true;
            opts.tfce.stat = 'tfce';
            a = a + 1;
            
        case '-Tstat', % not in the help
            
            % Type of cluster statistic
            opts.tfce.stat = vararginx{a+1};
            if ~ any(strcmp(opts.tfce.stat,{'tfce','density','tippett'})),
                error('TFCE statistic "%s" unknown.',opts.tfce.stat);
            end
            a = a + 2;
            
        case '-Tuni', % advanced
            
            % Do TFCE for uni?
            opts.tfce.uni.do = true;
            a = a + 1;
            
        case '-Tnpc', % advanced
            
            % Do TFCE for NPC?
            opts.NPC = true;
            opts.tfce.npc.do = true;
            a = a + 1;
            
        case '-Tmv', % advanced
            
            % Do TFCE for MV?
            opts.MV = true;
            opts.tfce.mv.do = true;
            a = a + 1;
            
        case {'-tfce1D','-tfce1d'}, % basic
            
            % Shortcut for -tfce_H 2 -tfce_E 2 -tfce_C 6,
            % i.e., parameters for TFCE in 2D mode
            opts.tfce.H      = 2;
            opts.tfce.E      = 2;
            opts.tfce.conn   = 6;
            a = a + 1;
            
        case {'-tfce2D','-tfce2d'}, % basic
            
            % Shortcut for -tfce_H 2 -tfce_E 1 -tfce_C 26,
            % i.e., parameters for TFCE in 2D mode
            opts.tfce.H      = 2;
            opts.tfce.E      = 1;
            opts.tfce.conn   = 26;
            a = a + 1;
            
        case {'-tfce_H','-tfce_h'}, % advanced
            
            % TFCE H parameter
            opts.tfce.H = vararginx{a+1};
            if ischar(opts.tfce.H),
                opts.tfce.H = str2double(opts.tfce.H);
            end
            a = a + 2;
            
        case {'-tfce_E','-tfce_e'}, % advanced
            
            % TFCE E parameter
            opts.tfce.E = vararginx{a+1};
            if ischar(opts.tfce.E),
                opts.tfce.E = str2double(opts.tfce.E);
            end
            a = a + 2;
            
        case {'-tfce_C','-tfce_c'}, % advanced
            
            % TFCE connectivity
            opts.tfce.conn = vararginx{a+1};
            if ischar(opts.tfce.conn),
                opts.tfce.conn = str2double(opts.tfce.conn);
            end
            a = a + 2;
            
        case '-tfce_dh', % advanced
            
            % TFCE delta-h parameter
            opts.tfce.deltah = vararginx{a+1};
            if ischar(opts.tfce.deltah),
                if strcmpi(opts.tfce.deltah,'auto'),
                    opts.tfce.deltah = 0;
                else
                    opts.tfce.deltah = str2double(opts.tfce.deltah);
                end
            end
            a = a + 2;
            
        case '-tableasvolume', % basic
            
            % Treat tables (e.g., CSV inputs) as volume, such that TFCE can
            % be calculated. This is useful for TFCE over timeseries.
            opts.tableasvolume = true;
            a = a + 1;
            
        case '-within', % basic
            
            % Define whether should permute blocks as a whole or not
            opts.within = true;
            a = a + 1;
            
        case '-whole', % basic
            
            % Define whether should permute blocks as a whole or not
            opts.whole = true;
            a = a + 1;
            
        case '-ee', % basic
            
            % Exchangeable errors (EE)?
            % If yes, this means permutations.
            opts.EE = true;
            a = a + 1;
            
        case '-ise', % basic
            
            % Independent and symmetric errors (ISE)?
            % If yes, this means sign-flippings.
            opts.ISE = true;
            a = a + 1;
            
        case '-cmcp', % advanced
            
            % Define whether Conditional Monte Carlo should be used or not,
            % that is, ignoring repeated elements in the permutation set.
            opts.cmcp = true;
            a = a + 1;
            
        case '-cmcx', % advanced
            
            % Define whether repeated rows in X should be ignored or not
            % when defining the permutations, which constitutes another
            % form of CMC
            opts.cmcx = true;
            a = a + 1;
            
        case '-twotail', % basic
            
            % Do a two-tailed test for all t-contrasts?
            opts.twotail = true;
            a = a + 1;
            
        case '-concordant', % basic
            
            % For the NPC, favour alternatives with the same sign?
            opts.concordant = true;
            a = a + 1;
            
        case '-reversemasks', % basic
            
            % Reverse masks.
            opts.reversemasks = true;
            a = a + 1;
            
        case '-corrmod', % basic
            
            % Correct over modalities.
            opts.corrmod = true;
            a = a + 1;
            
        case '-corrcon', % basic
            
            % Correct over contrasts.
            opts.corrcon = true;
            a = a + 1;
            
        case '-saveparametric', % advanced
            
            % If the user wants to have also the parametric p-values.
            opts.savepara = true;
            a = a + 1;
            
        case '-saveglm', % advanced
            
            % If the user wants, save COPE and VARCOPEs in the 1st
            % permutation.
            opts.saveglm = true;
            a = a + 1;
            
        case {'-savecdf','-save1-p'}, % basic
            
            % Save 1-p values (CDF) instead of the P-values
            opts.savecdf = true;
            a = a + 1;
            
        case '-logp',  % basic
            
            % Convert the P-values or (1-P)-values to -log10 before saving.
            opts.savelogp = true;
            a = a + 1;
            
        case '-savemask', % advanced
            
            % If the user wants to have also the masks used for each.
            % modality
            opts.savemask = true;
            a = a + 1;
            
        case '-rmethod', % advanced
            
            % Which method to use for the regression/permutation?
            if narginx > a,
                methlist = {           ...
                    'Draper-Stoneman', ...
                    'Still-White',     ...
                    'Freedman-Lane',   ...
                    'terBraak',        ...
                    'Kennedy',         ... % Kennedy won't be in the help
                    'Manly',           ...
                    'Huh-Jhun',        ...
                    'Dekker'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Regression/Permutation method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.rmethod = methlist{methidx};
            else
                error([...
                    'The option -rmethod requires a method to be specified.\n'...
                    'Consult the documentation.%s'],'');
            end
            
        case '-npc' % basic
            
            % This is a shortcut to enable NPC with the default settings.
            opts.NPC = true;
            opts.npcmod = true;
            a = a + 1;
            
        case '-npcmethod', % basic
            
            % Do the non-parametric combination?
            if nargin == a || (nargin > a && strcmp(vararginx{a+1}(1),'-')),
                error('The option "-npcmethod" requires a combining method to be indicated.');
                
            elseif nargin > a,
                
                % Which combining function to use for the combination?
                methlist = {               ...
                    'Tippett',             ...
                    'Fisher',              ...
                    'Pearson-David',       ...
                    'Stouffer',            ...
                    'Wilkinson',           ...
                    'Winer',               ...
                    'Edgington',           ...
                    'Mudholkar-George',    ...
                    'Friston',             ...
                    'Darlington-Hayes',    ...
                    'Zaykin',              ...
                    'Dudbridge-Koeleman',  ...
                    'Dudbridge-Koeleman2', ...
                    'Taylor-Tibshirani',   ...
                    'Jiang'};
                methidx = strcmpi(vararginx{a+1},methlist);
                
                % Check if method exists, and load extra parameters if needed
                if ~any(methidx);
                    error('Combining method "%s" unknown.',vararginx{a+1});
                elseif any(strcmpi(vararginx{a+1},{...
                        'Wilkinson',       ...
                        'Zaykin',          ...
                        'Jiang'})),
                    if ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm = 0.05;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        a = a + 3;
                        plm.npcparm = eval(vararginx{a+2});
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif any(strcmpi(vararginx{a+1},{...
                        'Darlington-Hayes',   ...
                        'Dudbridge-Koeleman', ...
                        'Jiang'})),
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        plm.npcparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif strcmpi(vararginx{a+1},'Friston'),
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        plm.npcparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif strcmpi(vararginx{a+1},'Dudbridge-Koeleman2'),
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm  = 1;
                        plm.npcparm2 = 0.05;
                        a = a + 2;
                    else
                        if ischar(vararginx{a+2}),
                            plm.npcparm = eval(vararginx{a+2});
                        else
                            plm.npcparm = vararginx{a+2};
                        end
                        if nargin == a + 2 || ...
                                ischar(vararginx{a+3}) && ...
                                strcmpi(vararginx{a+3}(1),'-'),
                            plm.npcparm2 = 0.05;
                        elseif ischar(vararginx{a+3}),
                            plm.npcparm2 = eval(vararginx{a+3});
                        else
                            plm.npcparm2 = vararginx{a+3};
                        end
                        a = a + 4;
                    end
                else
                    a = a + 2;
                end
                opts.npcmethod = methlist{methidx};
            end
            
        case '-npcmod', % basic
            
            % NPC over modalities.
            opts.NPC    = true;
            opts.npcmod = true;
            a = a + 1;
            
        case '-npccon', % basic
            
            % NPC over contrasts -- that is, all contrasts, even contrasts
            % in different designs (if more than one -d is supplied).
            opts.NPC       = true;
            opts.npccon    = true;
            opts.syncperms = true;
            a = a + 1;
            
        case '-mv', % basic
            
            % Compute classic multivariate statistics
            if nargin == a,
                opts.MV = true;
                a = a + 1;
                
            elseif nargin > a && strcmp(vararginx{a+1}(1),'-'),
                opts.MV = true;
                a = a + 1;
                
            elseif nargin > a,
                
                % Which multivariate statistic to use?
                methlist = {            ...
                    'auto',             ...
                    'HotellingTsq',     ...
                    'Wilks',            ...
                    'Lawley',           ...
                    'Lawley-Hotelling', ...
                    'Pillai',           ...
                    'Roy',              ...
                    'Roy-ii',           ...
                    'Roy-iii',          ...
                    'CCA',              ...
                    'PLS'};
                methidx = strcmpi(vararginx{a+1},methlist);
                
                % Check if method exists, and load extra parameters if needed
                if ~any(methidx);
                    error('Multivariate statistic "%s" unknown.',vararginx{a+1});
                elseif strcmpi(vararginx{a+1},'CCA'),
                    opts.MV  = false;
                    opts.CCA = true;
                    opts.PLS = false;
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        opts.ccaorplsparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        opts.ccaorplsparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        opts.ccaorplsparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif strcmpi(vararginx{a+1},'PLS'),
                    opts.MV  = false;
                    opts.CCA = false;
                    opts.PLS = true;
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        opts.ccaorplsparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        opts.ccaorplsparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        opts.ccaorplsparm = vararginx{a+2};
                        a = a + 3;
                    end
                else
                    opts.MV  = true;
                    opts.CCA = false;
                    opts.PLS = false;
                    a = a + 2;
                end
                opts.mvstat = methlist{methidx};
            end
            
        case '-fdr', % basic
            
            % Compute FDR-adjusted p-values
            opts.FDR = true;
            a = a + 1;
            
        case {'-accel','-approx'}, % advanced
            
            % Choose a method to do the approximation of p-values
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-'),
                methlist = {   ...
                    'negbin',  ...
                    'tail',    ...
                    'noperm',  ...
                    'gamma',   ...
                    'lowrank'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~ any(methidx);
                    error('Approximation method "%s" unknown.',vararginx{a+1});
                end
                for mm = 1:numel(methlist),
                    opts.accel.(methlist{mm}) = methidx(mm);
                end
                
                % Extra parameters
                if opts.accel.negbin,
                    
                    % Number of exceedances:
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-'),
                        if ischar(vararginx{a+2}),
                            opts.accel.negbin = str2double(vararginx{a+2});
                        else
                            opts.accel.negbin = vararginx{a+2};
                        end
                        a = a + 3;
                    else
                        opts.accel.negbin = opts.accel.negbin_nexced;
                        a = a + 2;
                    end
                    
                elseif opts.accel.tail,
                    
                    % Define whether include or not the unpermuted stat:
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-'),
                        if ischar(vararginx{a+2}),
                            if     any(strcmpi(vararginx{a+2},{'out','G1out','T1out','true', '1'})),
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false; % defensive, as the uncorrected will be invalid here.
                            elseif any(strcmpi(vararginx{a+2},{'in', 'G1in', 'T1in', 'false','0'})),
                                opts.accel.G1out = false;
                            end
                        else
                            if vararginx{a+2},
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false; % defensive, as the uncorrected will be invalid here.
                            else
                                opts.accel.G1out = false;
                            end
                        end
                        a = a + 3;
                    else
                        a = a + 2;
                    end
                    
                elseif opts.accel.gamma,
                    
                    % Define whether include or not the unpermuted stat:
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-'),
                        if ischar(vararginx{a+2}),
                            if     any(strcmpi(vararginx{a+2},{'out','G1out','T1out','true', '1'})),
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false;
                            elseif any(strcmpi(vararginx{a+2},{'in', 'G1in', 'T1in', 'false','0'})),
                                opts.accel.G1out = false; % defensive, as the uncorrected will be invalid here.
                            end
                        else
                            if vararginx{a+2},
                                opts.accel.G1out = true;
                                opts.saveuncorrected = false; % defensive, as the uncorrected will be invalid here.
                            else
                                opts.accel.G1out = false;
                            end
                        end
                        a = a + 3;
                    else
                        a = a + 2;
                    end
                    
                elseif opts.accel.lowrank,
                    
                    % Fraction of voxels to be sampled (if < 1) or actual
                    % number of voxels to be sampled.
                    if narginx > a+1 && ~strcmpi(vararginx{a+2}(1),'-'),
                        if ischar(vararginx{a+2}),
                            opts.accel.lowrank_val = str2double(vararginx{a+2});
                        else
                            opts.accel.lowrank_val = vararginx{a+2};
                        end
                        a = a + 3;
                    else
                        a = a + 2;
                    end
                    
                else
                    a = a + 2;
                end
            else
                error([...
                    'The options "-accel" and "-approx" require a method to.\n' ...
                    'be specified. Consult the documentation.%s'],'');
            end
            
        case {'-noniiclass','-nonifticlass'} % advanced
            
            % Disable using the NIFTI class
            opts.useniiclass = false;
            a = a + 1;
            
        case '-precision', % advanced
            
            % Precision to use?
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-'),
                methlist = {'single','double'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Precision "%s" unknown. Use "single" or "double".',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.precision = methlist{methidx};
            else
                error([...
                    'The option "-precision" requires a method to be specified.\n'...
                    'Use "-precision double" or "-precision single".']);
            end
            
        case '-saveperms', % advanced
            
            % Save the permutations
            opts.saveperms = true;
            a = a + 1;
            
        case '-savemetrics', % advanced
            
            % Save a file with the number of permutations, average
            % Hamming distance, etc.
            opts.savemetrics = true;
            a = a + 1;
            
        case '-inormal', % advanced
            
            % Inverse-normal transformation?
            opts.inormal = true;
            
            % Take the parameters given to -inormal
            parms = {};
            if narginx - a >= 1,
                if ~strcmp(vararginx{a+1}(1),'-'),
                    parms{1} = vararginx{a+1};
                end
            end
            if narginx - a >= 2,
                if ~strcmp(vararginx{a+2}(1),'-'),
                    parms{2} = vararginx{a+2};
                end
            end
            a = a + 1 + numel(parms);
            
            % Then modify the variables accordingly
            methlist = {   ...
                'Blom',    ...
                'Tukey',   ...
                'Bliss',   ...
                'Waerden', ...
                'SOLAR'};
            for p = 1:numel(parms),
                methidx = strcmpi(parms{p},methlist);
                if any(methidx);
                    opts.inormal_meth = parms{p};
                elseif any(parms{p},{'quali','qualitative','discrete'}),
                    opts.inormal_quanti = false;
                elseif any(parms{p},{'quanti','quantitative','continuous'}),
                    opts.inormal_quanti = true;
                else
                    error('Parameter "%s" unknown for the "-inormal" option.',parms{p});
                end
            end
            
        case '-probit', % advanced
            
            % Probit transformation?
            opts.probit = true;
            a = a + 1;
            
        case '-seed', % advanced
            
            % Seed for the random number generator
            opts.seed = vararginx{a+1};
            if ischar(opts.seed) && ...
                    ~any(strcmpi(opts.seed,{'shuffle','twist','reset'})),
                opts.seed = str2double(opts.seed);
            end
            a = a + 2;
            
        case '-demean', % basic
            
            % Demean data and design. Additionally, remove
            % a global intercept, if any, from the design.
            opts.demean = true;
            a = a + 1;
            
        case '-vgdemean', % advanced
            
            % Demean data and design within VG. Additionally, remove
            % a global intercept, if any, from the design.
            opts.vgdemean = true;
            a = a + 1;
            
        case '-ev4vg', % advanced
            
            % Add to the design matrix one EV for each variance group.
            opts.ev4vg = true;
            a = a + 1;
            
        case '-removevgbysize', % advanced
            
            % Remove from the analysis observations that are the only
            % in their variance group.
            opts.removevgbysize = vararginx{a+1};
            if ischar(opts.removevgbysize),
                opts.removevgbysize = str2double(opts.removevgbysize);
            end
            a = a + 2;
            
        case '-zstat', % advanced
            
            % Convert the statistic for each test to a z-score
            opts.zstat = true;
            a = a + 1;
            
        case '-pearson', % basic
            
            % Compute the Pearson's correlation coefficient (R^2 if rank(C)>1).
            opts.pearson = true;
            a = a + 1;
            
        case '-noranktest', % advanced
            
            % Don't test the rank(Y) before doing MANOVA/MANCOVA.
            opts.noranktest = true;
            a = a + 1;
            
        case '-transposedata', % advanced
            
            % Transpose the data if it's 2D?
            opts.transposedata = true;
            a = a + 1;
            
        case '-inputmv', % advanced
            
            % Treat the (sole) input as multivariate, that is,
            % each column is a variable in a multivariate model,
            % as opposed to independent univariate tests.
            % Useful with non-imaging data.
            opts.inputmv = true;
            a = a + 1;
            
        case '-verbosefilenames', % advanced
            
            % Don't simplify filenames when otherwise it would be possible
            opts.verbosefilenames = true;
            a = a + 1;
            
        case '-syncperms', % advanced
            
            % Synchronize permutations regardless of other options?
            opts.syncperms = true;
            a = a + 1;
            
        case {'-designperinput','-erie'}, % advanced
            
            % Use one design matrix for each input dataset?
            % This enables
            % syncP regardless.
            opts.designperinput = true;
            opts.syncperms      = true;
            a = a + 1;
            
        case '-subjidx', % advanced
            
            % Indices of the subjects to keep in the design
            opts.subjidx = vararginx{a+1};
            a = a + 2;
            
        case '-quiet', % basic
            
            % Don't print lines indicating progress
            opts.showprogress = false;
            a = a + 1;
            
        case '-nounivariate', % advanced
            
            % Save or not univariate tests? Useful with MV/NPC/CCA
            opts.saveunivariate = false;
            a = a + 1;
            
        case '-nouncorrected',  % advanced
            
            % Save or not uncorrected p-vals
            opts.saveuncorrected = false;
            a = a + 1;
            
        case '-saveuncorrected',  % advanced
            
            % Save or not uncorrected p-vals
            opts.saveuncorrected = true;
            a = a + 1;
            
        case '-savedof', % advanced
            
            % Save a file with the degrees of freedom?
            opts.savedof = true;
            a = a + 1;
            
        case '-pmethodp', % advanced
            
            % Which method to use to partition the model when defining
            % the permutations?
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-'),
                methlist = {    ...
                    'Guttman',  ...
                    'Beckmann', ...
                    'Ridgway',  ...
                    'none'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Partition method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.pmethodp = methlist{methidx};
            else
                error([...
                    'The option "-pmethodp" requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end
            
        case '-pmethodr', % advanced
            
            % Which method to use to partition the model when defining
            % doing the actual regression?
            if narginx > a && ~strcmpi(vararginx{a+1}(1),'-'),
                methlist = {    ...
                    'Guttman',  ...
                    'Beckmann', ...
                    'Ridgway',  ...
                    'none'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Partition method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.pmethodr = methlist{methidx};
            else
                error([...
                    'The option "-pmethodr" requires a method to be specified.\n' ...
                    'Consult the documentation.']);
            end
            
        case '-forceintersectionmask', % currently not in the help
            
            % Force the generation of an intersection mask, even if there
            % is no MV or NPC (useful for mediation analysis).
            opts.forcemaskinter = true;
            a = a + 1;
            
        otherwise
            error('Unknown option: "%s"',vararginx{a});
    end
end

% Check for the existence of other programs for input/output
palm_checkprogs;

if Ni == 0,
    error('Missing input data (missing "-i").');
end

% Make sure the NPC options make sense
% - if -npc only, it's -npcmod
% - if -npcmod only, it's -npcmod
% - if -npccon only, it's -npccon
% - if -npcmod and -npccon, it's both
if opts.NPC && ~ opts.npccon,
    opts.npcmod = true;
end

% A quick check for the case of 1 EV per column of Y.
if opts.evperdat,
    if any([...
            opts.ev4vg
            opts.pearson]),
        error([...
            'The option "-evperdat" is incompatible with the options listed below:\n' ...
            '"-ev4vg"\n' ...
            '"-pearson"%s'],'');
    end
    if strcmpi(opts.rmethod,'terBraak'),
        error('The option "-evperdat" cannot be used with the ter Braak method (not implemented)');
    end
end

% No spatial statistics if no univariate results will be saved anyway
if ~ opts.saveunivariate,
    opts.cluster.uni.do = false;
    opts.tfce.uni.do    = false;
end

% This simplifies some tests later
opts.spatial.do  = false;
opts.spatial.uni = false;
opts.spatial.npc = false;
opts.spatial.mv  = false;
if any([ ...
        opts.cluster.uni.do  ...
        opts.cluster.npc.do  ...
        opts.cluster.mv.do   ...
        opts.tfce.uni.do     ...
        opts.tfce.npc.do     ...
        opts.tfce.mv.do]'),
    opts.spatial.do = true;
    if any([ ...
            opts.cluster.uni.do  ...
            opts.tfce.uni.do]'),
        opts.spatial.uni = true;
    end
    if any([ ...
            opts.cluster.npc.do  ...
            opts.tfce.npc.do]'),
        opts.spatial.npc = true;
    end
    if any([ ...
            opts.cluster.mv.do   ...
            opts.tfce.mv.do]'),
        opts.spatial.mv = true;
    end
end

% Adjust the delta-h.
if any([ ...
        opts.tfce.uni.do
        opts.tfce.npc.do
        opts.tfce.mv.do ]) && ...
        strcmpi(opts.tfce.deltah,'auto'),
    opts.tfce.deltah = 0;
end

% Sanity checks for the acceleration modes.
if sum(logical([ ...
        opts.accel.negbin, ...
        opts.accel.tail,   ...
        opts.accel.noperm, ...
        opts.accel.gamma,  ...
        opts.accel.lowrank])) > 1,
    error('Only one approximation method can be used for a given run.');
end
if opts.accel.negbin,
    if opts.accel.negbin < 2,
        error('The parameter r given to "-accel negbin <r>" must be >= 2.')
    end
    if opts.nP0 < 3,
        error('The option "-accel negbin <r>" needs at least 3 permutations.')
    end
    if ~ opts.saveuncorrected,
        error('The option "-nouncorrected" cannot be used with "-accel negbin".');
    end
    if (opts.corrmod || opts.corrcon) && ~ opts.FDR,
        error('The option "-accel negbin" cannot be used with FWER-correction, only FDR.');
    end
    if opts.NPC,
        error('The option "-accel negbin" cannot be used with NPC.');
    end
    if opts.spatial.do,
        error('The option "-accel negbin" cannot be used with spatial statistics.');
    end
    if opts.saveperms,
        error('The option "-saveperms" cannot be used together with "-accel negbin".');
    end
end
if opts.accel.noperm,
    if ~ opts.saveuncorrected,
        error('The option "-nouncorrected" cannot be used with "-accel noperm".');
    end
    if (opts.corrmod || opts.corrcon) && ~ opts.FDR,
        error('The option "-accel noperm" cannot be used with FWER-correction, only FDR.');
    end
    if opts.NPC,
        error('The option "-accel noperm" cannot be used with NPC.');
    end
    if opts.spatial.do,
        error('The option "-accel noperm" cannot be used with spatial statistics.');
    end
    if opts.MV,
        if ~ any(strcmpi(opts.mvstat,{'auto','pillai'})),
            warning([...
                'With multivariate tests, the option "-accel noperm" can\n' ...
                '         only be used with the Pillai'' trace statistic.\n' ...
                '         Changing automatically to Pillai.%s'],'');
        end
        opts.mvstat = 'pillai';
    end
    if Ni > 1 && ~ opts.MV,
        error([...
            'The option "-accel noperm" needs to be used with a single modality\n'...
            '       modality or with "-mv".%s'],'');
    end
end
if opts.accel.lowrank,
    if opts.pearson,
        error('The option "-accel lowrank" cannot be used with "-pearson".');
    end
    if opts.MV,
        error('The option "-accel lowrank" cannot be used with MV.');
    end
    if opts.NPC,
        error('The option "-accel lowrank" cannot be used with NPC.');
    end
    if opts.spatial.do,
        error('The option "-accel lowrank" cannot be used with spatial statistics.');
    end
    if opts.evperdat,
        error('The option "-accel lowrank" cannot be used with "-evperdat".');
    end
    if opts.missingdata,
        error('The option "-accel lowrank" cannot be used with missing data.');
    end
        if opts.saveglm,
        error('The option "-accel lowrank" cannot be used with "-saveglm".');
    end
end

% Some options can't be used with missing data
if Nimiss || Ndmiss,
    opts.missingdata = true;
    if opts.MV,
        error('The option "-mv" cannot be used with missing data.');
    end
    if ~ opts.zstat && ~ opts.mcar,
        warning([...
            'With missing data MAR/MNAR, the option "-zstat" is mandatory.\n' ...
            '         Adding it automatically.%s'],'');
        opts.zstat = true;
    end
    if ~ opts.cmcx && ~ opts.mcar,
        warning([...
            'With missing data, the option "-cmcx" is mandatory.\n' ...
            '         Adding it automatically.%s'],'');
        opts.cmcx = true;
    end
    if opts.demean && ~ opts.mcar,
        warning([...
            'With missing data, the option "-demean" must not be used.\n' ...
            '         Removing it automatically.%s'],'');
        opts.demean = false;
    end
    if ~ strcmpi(opts.pmethodp,'guttman') || ~ strcmpi(opts.pmethodr,'guttman'),
        warning([...
            'With missing data, the partitioning must use the "Guttman".\n' ...
            '         method. Adding automatically the options\n' ...
            '         "-pmethodp Guttman" and "-pmethodr Guttman".%s'],'');
        opts.pmethodp = 'guttman';
        opts.pmethodr = 'guttman'; 
    end
    if opts.ev4vg || opts.removevgbysize > 0,
        error('Missing data cannot be used with "-ev4vg" or "-removevgbysize".')
    end
end

% Some NPC methods don't have an analytical form for the parametric p-value
if opts.NPC && any(strcmpi(opts.npcmethod,{'Darlington-Hayes','Jiang'})),
    plm.nonpcppara = true;
    if opts.savepara,
        warning([...
            'No parametric combination p-value will be saved for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
    end
    if opts.spatial.npc,
        warning([ ...
            'No NPC spatial statistic will be produced for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
        opts.cluster.npc.do = false;
        opts.tfce.npc.do    = false;
        opts.spatial.npc    = false;
    end
else
    plm.nonpcppara = false;
end

% Likewise, some MV methods don't have an analytical form for the parametric p-value
if opts.MV && strcmpi(opts.mvstat,'Roy_iii'),
    plm.nomvppara = true;
    if opts.savepara,
        warning('No parametric p-value will be saved for the Roy_iii method.%s','');
    end
    if opts.spatial.mv,
        warning([ ...
            'No multivariate cluster-level or TFCE statistic will be produced\n', ...
            '         for the Roy_iii statistic%s'],'');
        opts.cluster.mv.do = false;
        opts.tfce.mv.do    = false;
        opts.spatial.mv    = false;
    end
elseif (opts.CCA || opts.PLS) && opts.savepara,
    plm.nomvppara = true;
    warning([...
        'No parametric p-value will be saved for CCA or PLS.\n', ...
        '         Use Wilks'' lambda instead.%s'],'');
else
    plm.nomvppara = false;
end

% Some more warnings and sanity checks
if opts.savecdf && opts.savelogp,
    error('Should not use "-save1-p" together with "-logp"');
end
if ~ opts.inputmv && opts.designperinput && Ni ~= Nd,
    error([
        'To use the option "-designperinput", the number of design files must\n' ...
        'match the number of inputs.\n%s'],'');
end
if (Nt || Nf) && Ncon,
    error('Cannot mix options "-t" or "-f" with "-con".');
end
if Nt || Nf,
    if Nt > Nd,
        error('More t-contrast files (%d) than valid design files (%d) were supplied.',Nt,Nd);
    end
    if Nt ~= 1 && Nt ~= Nd,
        error(['The number of supplied t-contrast files (option "-t") must be 1 or\n',...
            'the same as the number of design files (option "-d") (%d).'],Nd);
    end
    if Nf > Nt,
        error('More F-contrast files (%d) than t-contrast files (%d) were supplied.',Nf,Nt);
    end
elseif Ncon,
    if Ncon > Nd,
        error('More contrast files (%d) than design files (%d) were supplied.',Nt,Nd);
    end
    if Ncon ~= 1 && Ncon ~= Nd,
        error(['The number of supplied contrast files (option "-con") must be 1 or\n',...
            'the same as the number of design files (option "-d") (%d).'],Nd);
    end
end
if opts.pearson && (opts.NPC || opts.MV),
    error([ ...
        'It isn''t possible to compute the Pearson''s r or R^2 together with NPC or\n', ...
        '         multivariate methods.%s'],'');
end
if opts.pearson && ~ opts.demean,
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the data and the design\n' ...
        '         must be mean centered. Adding option "-demean".%s'],'');
    opts.demean = true;
end
if opts.pearson && opts.savepara,
    error([ ...
        'The option "-saveparametric" cannot be used with "-pearson".\n' ...
        '         For a parametric p-value, drop "-pearson" and use the\n' ...
        '         respective results for the t or F-statistics.%s'],'');
end
if opts.pearson && ~ any(strcmpi(opts.pmethodr,{'beckmann','ridgway'})),
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the design must be\n' ...
        '         partitioned using the Beckmann or Ridgway schemes.'...
        '         Adding the option "-pmethodr Beckmann".%s'],'');
    opts.pmethodr = 'beckmann';
end
if (opts.CCA || opts.PLS) && ~ opts.demean,
    warning([ ...
        'To perform CCA or PLS, the data and the design\n' ...
        '         must be mean centered. Adding option "-demean".%s'],'');
    opts.demean = true;
end
if (opts.CCA || opts.PLS) && ~ any(strcmpi(opts.pmethodr,{'beckmann','ridgway'})),
    warning([ ...
        'To perform CCA or PLS, the design must be\n' ...
        '         partitioned using the Beckmann or Ridgway schemes.'...
        '         Adding the option "-pmethodr Beckmann".%s'],'');
    opts.pmethodr = 'beckmann';
end
if opts.demean && opts.vgdemean && ~opts.pearson && ~opts.CCA && ~opts.PLS,
    warning([...
        'Cannot use the option "-demean" together with "-vgdemean"\n'...
        '         Ignoring the option "-vgdemean".%s'],'');
    opts.vgdemean = false;
end
if opts.ev4vg && opts.vgdemean,
    error('Cannot use the option "-ev4vg" together with "-vgdemean"');
end
if opts.MV && ~ opts.singlevg,
    error('Currently MV is only allowed with a single VG, that is, assuming homoscedasticity.');
end
if opts.designperinput && opts.MV,
    error('It''s not possible to use the option "-designperinput" together with the option "-mv".');
end
if ~opts.saveunivariate && ~opts.MV && ~opts.NPC && ~opts.CCA && ~opts.PLS,
    error('The option "-nounivariate" can only be used with "-mv" or "-npc".');
end
if opts.MV && (opts.CCA || opts.PLS),
    error('Cannot do classical MANOVA at the same time as CCA or PLS.');
end
if Ni > 1 && opts.inputmv,
    error('Option "-inputmv" cannot be used with more than one "-i".')
end
if opts.concordant && ~ opts.NPC,
    error('The option "-concordant" is for use with NPC only.');
end
if opts.concordant && opts.twotail,
    error(['Cannot use "-concordant" together with "-twotail" (inadmissible).\n'...
        'Use either of these, but not both together.%s'],'');
end
if opts.tonly && opts.fonly,
    error('Cannot use "-tonly" together with "-fonly".');
end
if opts.probit && opts.inormal,
    error('Cannot use "-probit" together with "-inormal".');
end
if opts.FDR && ~ opts.saveuncorrected,
    error(['Option "-fdr" cannot be used together with "-nouncorrected".\n'...
        'Use either of these, but not both together. If you are using tail or\n'...
        'gamma approximations, consider keeping "-nouncorrected", and leave\n'...
        '"-fdr" for a separate call in which tail or gamma are not used.%s'],'');
end

% Initialize the random number generator (if nP = 0, no need for that)
if opts.nP0,
    if palm_isoctave,
        if any(strcmpi(opts.seed,{'reset','shuffle','twist'})),
            opts.seed = 'reset';
        end
        rand('state',opts.seed); %#ok
    else
        if any(strcmpi(opts.seed,{'reset','shuffle','twist'})),
            opts.seed = 'shuffle';
        end
        rng('default');
        rng(opts.seed);
    end
end

% Read and organise the surfaces. If no surfaces have been loaded, but the
% user wants cluster extent, cluster mass, or TFCE, an error will be
% printed later down in the code.
% At this stage the average area from native geometry is also loaded if it
% was provided. If a weight was given, such as 1, use this weight, so that
% all faces are treated as if having the same size. If nothing was
% supplied, use the actual area of the surface.
if opts.spatial.do && Ns > 0,
    plm.srf     = cell(Ns,1);
    plm.srfarea = cell(Ns,1);
    plm.srfadj  = cell(Ns,1);
    for s = 1:Ns,
        
        % Load surface
        fprintf('Loading surface %d/%d: %s\n',s,Ns,opts.s{s});
        plm.srf{s} = palm_miscread(opts.s{s});
        
        % Load areas
        if isempty(opts.sa{s}),
            % No area means area of the actual surface file
            plm.srfarea{s}.data = [];
        elseif exist(opts.sa{s},'file'),
            % A file with the average areas from native geometry
            plm.srfarea{s} = palm_miscread(opts.sa{s},opts.useniiclass,opts.o,opts.precision);
        elseif ~ isnan(str2double(opts.sa{s})),
            % A weight (such as 1)
            plm.srfarea{s}.data = str2double(opts.sa{s});
        else
            % Otherwise gives a helpful error message:
            error('Unknown option given to "-s" or file doesn''t exist:\n%s',opts.sa{s});
        end
    end
end

% There should be no more masks than modalities, and the number of
% masks needs to be either 1 or the same number of modalities.
if Nm > Ni,
    error([...
        'There are more masks supplied with "-m" (%d masks) than\n'...
        'modalities supplied with "-i" (%d modalities)'],Nm,Ni);
elseif Nm > 1 && Nm ~= Ni,
    error([...
        'The number of masks supplied with "-m" (%d masks) is larger than 1,\n'...
        'but still not the same as the number of modalities supplied with\n'...
        'the option "-i" (%d modalities).'],Nm,Ni);
end

% Read and organise the masks. If there are no masks specified, one for
% each modality will be created after each modality is loaded.
plm.masks = cell(Ni,1);
for m = 1:Nm,
    plm.masks{m} = palm_miscread(opts.m{m},opts.useniiclass,opts.o,opts.precision);
    if strcmp(plm.masks{m}.readwith,'nifticlass'),
        plm.masks{m}.data = double(plm.masks{m}.data);
    end
    if opts.reversemasks,
        plm.masks{m}.data(isnan(plm.masks{m}.data)) = 1;
        plm.masks{m}.data(isinf(plm.masks{m}.data)) = 1;
        plm.masks{m}.data = ~ logical(plm.masks{m}.data);
    else
        plm.masks{m}.data(isnan(plm.masks{m}.data)) = 0;
        plm.masks{m}.data(isinf(plm.masks{m}.data)) = 0;
        plm.masks{m}.data = logical(plm.masks{m}.data);
    end
end
if Nm == 1,
    for i = 2:Ni,
        plm.masks{i} = plm.masks{1};
    end
end

% Indices of the subjects that will be kept
if ~ isempty(opts.subjidx),
    plm.subjidx = palm_miscread(opts.subjidx);
    plm.subjidx = round(plm.subjidx.data);
end

% Read and organise the data.
plm.Yset     = cell(Ni,1);  % Regressands (Y)
plm.Yisvol   = false(Ni,1); % Is Y a volume image?
plm.Yissrf   = false(Ni,1); % Is Y a surface-based image (DPX)?
plm.Yisvtx   = false(Ns,1); % Is vertexwise?
plm.Yisfac   = false(Ns,1); % is facewise? (this is currently dichotomous with Yisvtx, but later there may be edges/connectivity too)
plm.Yarea    = cell(Ns,1);  % To store area per face or per vertex (used for cluster-level & TFCE).
plm.Ykindstr = cell(Ni,1);  % string to save the files later
for i = 1:Ni,
    
    % Read input file
    fprintf('Reading input %d/%d: %s\n',i,Ni,opts.i{i});
    if Nm == 0,
        maskstruct = [];
    else
        maskstruct = plm.masks{i};
    end
    [plm.Yset{i},plm.masks{i},plm.Yisvol(i),plm.Yissrf(i),plm.Ykindstr{i}] = ...
        palm_ready(opts.i{i},maskstruct,opts,true);
    
    % Select subjects
    if ~ isempty(plm.subjidx),
        plm.Yset{i} = plm.Yset{i}(plm.subjidx,:);
    end
    
    % For the first input data, keep the size to
    % compare with the others, then check the size
    if i == 1,
        plm.N = size(plm.Yset{i},1);
    end
    if size(plm.Yset{i},1) ~= plm.N,
        error([
            'At least two of the input data files do not have\n' ...
            'compatible sizes:\n' ...
            '- File %d (%s) has %d observations\n' ...
            '- File %d (%s) has %d observations'], ...
            1,opts.i{1},plm.N, ...
            i,opts.i{i},size(plm.Yset{i},1));
    end
    
    % If this is a DPX/curvature file, and if one of the spatial
    % statistics has been invoked, check if surfaces are available
    % and with compatible size, then compute the area (dpv or dpf).
    % Also take this opportunity to compute the adjacency matrix.
    if opts.tableasvolume,
        plm.Yisvol(i) = true;
    else
        if opts.spatial.do && Ns > 0,
            
            
            % String defining the types, for the filenames and other tasks.
            if Ns == 1, s = 1; else s = i; end
            if any(size(plm.srf{s}.data.vtx,1) == ...
                    size(plm.masks{i}.data));
                plm.Yisvol(i)   = false;
                plm.Yissrf(i)   = true;
                plm.Yisvtx(i)   = true;
                plm.Yisfac(i)   = false;
                plm.Ykindstr{i} = '_dpv';
            elseif any(size(plm.srf{s}.data.fac,1) == ...
                    size(plm.masks{i}.data));
                plm.Yisvol(i)   = false;
                plm.Yissrf(i)   = true;
                plm.Yisvtx(i)   = false;
                plm.Yisfac(i)   = true;
                plm.Ykindstr{i} = '_dpf';
            else
                error([...
                    'Surface file does not match the input data:\n' ...
                    '- Surface file %d has %d vertices and %d faces (%s)\n' ...
                    '- Input data file %d has %d points (%s)'],...
                    s,size(plm.srf{s}.data.vtx,1),size(plm.srf{s}.data.fac,1),opts.s{s},...
                    i,max(size(plm.masks{i}.data)),opts.i{i});
            end
            
            % Surface area, to be used by the spatial statistics
            if isempty(plm.srfarea{s}.data),
                
                % No area file given, use the actual surface area
                plm.Yarea{i} = palm_calcarea(plm.srf{s}.data,plm.Yisvtx(i));
                
            elseif numel(plm.srfarea{s}.data) == 1,
                
                % A weight given (such as 1): use that for each vertex or face,
                % treating as if all had the same area.
                if plm.Yisvtx(i),
                    plm.Yarea{i} = plm.srfarea{s}.data .* ...
                        ones(size(plm.srf{s}.data.vtx,1),1);
                elseif plm.Yisfac(i),
                    plm.Yarea{i} = plm.srfarea{s}.data .* ...
                        ones(size(plm.srf{s}.data.fac,1),1);
                end
                
            else
                
                % Otherwise, just use the data from the file (already loaded).
                plm.Yarea{i} = plm.srfarea{s}.data;
            end
            
            % Compute the adjacency matrix
            plm.Yadjacency{i} = palm_adjacency(plm.srf{s}.data.fac,plm.Yisvtx(i));
            
        elseif opts.spatial.do && Ns == 0 && any(strcmpi(plm.masks{i}.readwith,{'load','fs_load_mgh','gifti'})),
            error([ ...
                'To use spatial statistics with vertexwise or facewise data it is\n'...
                'necessary to provide the surface files (with the option "-s").%s'],'');
        end
    end
end
plm.nmasks = numel(plm.masks);
plm.nY     = numel(plm.Yset); % this is redefined below if opts.inputmv is set.

% Some extra packages for Octave
if opts.spatial.do && palm_isoctave && any(plm.Yisvol),
    pkg load image
end

% Read the EV per datum:
if opts.evperdat,
    
    % Some sanity check:
    opts.evpos = cat(1,opts.evpos{:});
    if size(unique(opts.evpos,'rows'),1) ~= size(opts.evpos,1);
        error([
            'Some EV per datum have been defined for the same\n'...
            'position in the same design matrices.%s'],'');
    end
    plm.EVset  = cell(Nevd,1);
    plm.nEVdat = Nevd;
    
    % If there's one design per input, use the same masks as
    % those of the input files. Otherwise, create them.
    if ~ opts.designperinput && Nm == 1,
        for ev = 1:plm.nEVdat,
            plm.masksEV{ev} = plm.masks{1};
        end
    elseif opts.designperinput || (plm.nY == 1 && Nd == 1),
        for ev = 1:plm.nEVdat,
            plm.masksEV{ev} = plm.masks{opts.evpos(ev,2)};
        end
    else
        plm.masksEV = cell(plm.nEVdat,1);
    end
    
    % Read input file & select subjects
    for ev = 1:plm.nEVdat,
        fprintf('Reading EV per datum %d/%d: %s\n',ev,plm.nEVdat,opts.evdatfile{ev});
        [plm.EVset{ev},plm.masksEV{ev}] = palm_ready(opts.evdatfile{ev},plm.masksEV{ev},opts,false);
        if ~ isempty(plm.subjidx),
            plm.EVset{ev} = plm.EVset{ev}(plm.subjidx,:);
        end
    end
    plm.nmasksEV = numel(plm.masksEV);
    
    % Make the intersection of the EVs that will go all in the same design
    Dlist = unique(opts.evpos(:,2),'rows');
    for d = 1:numel(Dlist,1),
        evidx = find(opts.evpos(:,2) == Dlist(d));
        if numel(evidx) > 1,
            newmask = true(size(plm.masksEV{1}.data));
            for ev = evidx',
                newmask = newmask & plm.masksEV{ev}.data;
            end
            for ev = evidx',
                plm.EVset{ev} = plm.EVset{ev}(:,newmask(plm.masksEV{ev}.data(:)));
                plm.masksEV{ev}.data = newmask;
            end
        end
    end
    
    % Then make the intersections with the respective masks of the input files
    if (opts.designperinput || (plm.nY == 1 && Nd == 1)) && ...
            ~ opts.npcmod && ~ opts.npccon && ~ opts.MV,
        for d = unique(opts.evpos(:,2))',
            EV = find(opts.evpos(:,2) == d);
            newmask = zeros(numel(plm.masks{d}.data),numel(EV)+1);
            c = 1;
            for ev = EV',
                newmask(:,c) = plm.masksEV{ev}.data(:);
                c = c + 1;
            end
            newmask(:,c) = plm.masks{d}.data(:);
            newmask = reshape(all(newmask,2),size(plm.masks{d}.data));
            for ev = EV',
                plm.EVset{ev} = plm.EVset{ev}(:,newmask(plm.masksEV{ev}.data(:)));
                plm.masksEV{ev}.data = newmask;
            end
            plm.Yset{d} = plm.Yset{d}(:,newmask(plm.masks{d}.data(:)));
            plm.masks{d}.data    = newmask;
        end
    else
        sizy = zeros(plm.nmasks,1);
        for y = 1:plm.nmasks,
            sizy(y) = size(plm.masks{y},1);
        end
        sizev = zeros(plm.nmasksEV,1);
        for ev = 1:plm.nmasksEV,
            sizev(ev) = size(plm.masksEV{ev},1);
        end
        if  numel(unique(sizy)) > 1 || ...
                numel(unique(sizev)) > 1 || ...
                sizy(1) ~= sizev(1),
            error([...
                'For multiple "-i" and/or "-evperdat", with "-npccon",\n',...
                'and without the option "-designperinput", the inputs, EVs \n'...
                'and masks need to be all of the same sizes.%s'],'');
        end
        newmask = true(size(plm.masksEV{1}.data));
        for y = 1:plm.nmasks,
            newmask = newmask & plm.masks{y}.data;
        end
        for ev = 1:plm.nmasksEV,
            newmask = newmask & plm.masksEV{ev}.data;
        end
        for y = 1:plm.nmasks,
            plm.Yset{y} = plm.Yset{y}(:,newmask(plm.masks{y}.data(:)));
            plm.masks{y}.data = newmask;
        end
        for ev = 1:plm.nmasksEV,
            plm.EVset{ev} = plm.EVset{ev}(:,newmask(plm.masksEV{ev}.data(:)));
            plm.masksEV{ev}.data = newmask;
        end
    end
end
clear('newmask');

% Create an intersection mask if NPC or MV is to be done, and further apply
% to the data that was previously masked above, as needed.
if opts.npcmod || opts.MV || opts.CCA || opts.PLS || opts.forcemaskinter,
    if plm.nmasks > 1,
        
        % If there is one mask per modality, make an instersection mask.
        maskinter = true(size(plm.masks{1}.data));
        for m = 1:plm.nmasks,
            maskinter = maskinter & plm.masks{m}.data;
        end
        if opts.evperdat,
            for ev = 1:plm.nEVdat,
                maskinter = maskinter & plm.masksEV{ev}.data;
            end
        end
        
        % Note that this line below uses Ytmp, which is from the previous loop.
        % This can be used here because with NPC all data has the same size.
        plm.maskinter = palm_maskstruct(maskinter(:)',plm.masks{1}.readwith,plm.masks{1}.extra);
        
        % Apply it to further subselect data points
        for y = 1:plm.nY,
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.masks{y}.data));
        end
        if opts.evperdat,
            for ev = 1:plm.nEVdat,
                plm.EVset{ev} = plm.EVset{ev}(:,plm.maskinter.data(plm.masksEV{ev}.data));
            end
        end
    else
        
        % If only one mask was given.
        plm.maskinter = plm.masks{1};
        for y = 1:plm.nY,
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.masks{1}.data));
        end
    end
end

% Sizes for later
if opts.evperdat,
    plm.EVsiz = zeros(plm.nEVdat,1);
    for ev = 1:plm.nEVdat,
        plm.EVsiz(ev) = size(plm.EVset{ev},2);
    end
end

% Make sure that all data have the same size if NPC or MV are to be done
if opts.npcmod || opts.MV || opts.forcemaskinter,
    % The plm.Ysiz is redefined below.
    for y = 1:plm.nY,
        plm.Ysiz(y) = size(plm.Yset{y},2);
    end
    [usiz,uidx] = unique(plm.Ysiz);
    if numel(usiz) > 2 || (numel(usiz) == 2 && min(usiz) ~= 1),
        error('The sizes of some of the imaging modalities don''t match');
    elseif numel(usiz) == 2 && usiz(1) == 1,
        for y = 1:plm.nY,
            if plm.Ysiz(y) == 1,
                fprintf('Expanding modality #%d to match the size of the others.\n',y);
                plm.Yset{y} = repmat(plm.Yset{y},[1 usiz(2)]);
                if plm.nmasks > 1,
                    if numel(plm.masks{y}.data) == 1,
                        plm.masks{y}.data = plm.masks{uidx(2)}.data;
                    else
                        error('Modality expansion is only allowed for single input variables.')
                    end
                end
            end
        end
    end
    clear('usiz');
end

% Make sure none of the modalities is empty
for y = 1:plm.nY,
    if any(size(plm.Yset{y}) == 0),
        error('Modality %d has no data.\n',y);
    end
end

% If the multiple columns of the (sole) input are to be treated
% in a multivariate fashion
if opts.inputmv,
    tmp1 = plm.Yset{1};
    tmp2 = plm.Ykindstr{1};
    nY   = size(tmp1,2);
    plm.Yset     = cell(nY,1);
    plm.Ykindstr = cell(nY,1);
    for y = 1:nY,
        plm.Yset{y}     = tmp1(:,y);
        plm.Ykindstr{y} = tmp2;
    end
    clear tmp1 tmp2 nY;
    plm.nY = numel(plm.Yset);
end

% A variable with the cumulative sizes of all modalities will be handy later
plm.Ysiz = zeros(plm.nY,1);
for y = 1:plm.nY,
    plm.Ysiz(y) = size(plm.Yset{y},2);
end
plm.Ycumsiz = vertcat(0,cumsum(plm.Ysiz));

% Take this opportunity to save the masks if the user requested.
if opts.savemask,
    for y = 1:plm.nmasks,
        M = plm.masks{y};
        if plm.nY == 1 || plm.nmasks == 1,
            M.filename = sprintf('%smask',opts.o);
        else
            M.filename = sprintf('%smask_m%d',opts.o,y);
        end
        M.data = double(M.data);
        palm_miscwrite(M,true);
    end
    if plm.nY > 1 && (opts.npcmod || opts.MV || opts.forcemaskinter),
        M          = plm.maskinter;
        M.filename = sprintf('%sintersection_mask',opts.o);
        M.data     = double(M.data);
        palm_miscwrite(M,true);
    end
end

% If MV was selected, make sure that Y is full rank.
if opts.MV && ~ opts.noranktest,
    fprintf('Testing rank of the data for the MV tests. To skip, use -noranktest.\n')
    Y = cat(3,plm.Yset{:});
    Y = permute(Y,[1 3 2]);
    failed = false(1,size(Y,3));
    for v = 1:size(Y,3),
        if rank(Y(:,:,v)) ~= plm.nY;
            failed(v) = true;
        end
    end
    if any(failed),
        fname = sprintf('%smv_illconditioned',opts.o);
        palm_quicksave(double(failed),0,opts,plm,[],[],[],fname);
        error([
            'One or more datapoints have ill-conditioned data. It is\n' ...
            'not possible to run multivariate analyses as MANOVA/MANCOVA.\n' ...
            'Please, see these datapoints marked as 1 in the file:\n' ...
            '%s.*\n'],fname); %#ok
    end
end; clear Y;

% Applies an inverse-normal transformation to the modalities if the user requested
if opts.inormal,
    for y = 1:plm.nY,
        plm.Yset{y} = palm_inormal( ...
            plm.Yset{y},            ...
            opts.inormal_meth,      ...
            opts.inormal_quanti);
    end
end

% Applies a probit transformation to the modalities if the user requested
if opts.probit,
    for y = 1:plm.nY,
        if min(plm.Yset{y}(:)) < 0 || max(plm.Yset{y}(:)) > 1,
            error([
                'Probit transformation can only be used with data in the interval [0 1].\n' ...
                'This fails for at least modality #%d.'],y);
        end
        plm.Yset{y} = erfinv(2*(plm.Yset{y}*0.999999999999999 + .5e-15)-1)*sqrt(2);
    end
end

% Make the adjustments for the EE and ISE options.
% - if the user gives nothing, its EE by default.
% - if the user gives ISE only, it's ISE only
% - if the user gives EE only, it's EE only
% - if the user gives both, it's both
if ~opts.EE && ~opts.ISE,
    opts.EE  = true;
end

% Read and assemble the design matrices.
fprintf('Reading design matrix and contrasts.\n');
if opts.evperdat,
    plm.Mset = cell(max(Nd,max(opts.evpos(:,2))),1);
    for m = 1:numel(plm.Mset),
        plm.Mset{m} = ones(plm.N,1);
    end
else
    plm.Mset = cell(max(Nd,1),1);
end
plm.nM = numel(plm.Mset);
if Nd == 0 && ~ opts.evperdat,
    plm.Mset{1} = ones(plm.N,1);
    opts.EE     = false;
    opts.ISE    = true;
elseif Nd > 0,
    for m = 1:Nd,
        Mtmp = palm_miscread(opts.d{m},[],[],opts.precision);
        plm.Mset{m} = Mtmp.data;
        if ~ isempty(plm.subjidx) && size(plm.Mset{m},1) ~= plm.N,
            plm.Mset{m} = plm.Mset{m}(plm.subjidx,:);
        end
        if size(plm.Mset{m},1) ~= plm.N,
            error([
                'The number of rows in the design matrix does\n' ...
                'not match the number of observations in the data.\n' ...
                '- Rows in the matrix: %d\n' ...
                '- Observations in the data: %d\n' ...
                'In file %s\n'], ...
                size(plm.Mset{m},1),plm.N,opts.d{m});
        end
        if any(isnan(plm.Mset{m}(:))) || any(isinf(plm.Mset{m}(:))),
            error([
                'The design matrix cannot contain NaN or Inf.\n' ...
                'In file %s\n'],opts.d{m});
        end
    end
end

% Include the EV per datum
if opts.evperdat,
    for ev = 1:plm.nEVdat,
        if ndims(plm.Mset{opts.evpos(ev,2)}) == 2,
            plm.Mset{opts.evpos(ev,2)} = ...
                repmat(plm.Mset{opts.evpos(ev,2)},[1 1 plm.EVsiz(ev)]);
        end
        plm.Mset{opts.evpos(ev,2)}(:,opts.evpos(ev,1),:) = ...
            permute(plm.EVset{ev},[1 3 2]);
    end
    plm = rmfield(plm,{'EVset','EVsiz'});
end

% Some related sanity checks
if opts.evperdat,
    for m = 1:plm.nM,
        if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
        for y = loopY,
            if size(plm.Yset{y},2) ~= size(plm.Mset{m},3),
                error([
                    'The size of the data and the size of the EV per datum\n' ...
                    'don''t match.%s'],'')
            end
        end
    end
end

% Read and organise the contrasts for each design.
plm.Cset = cell(plm.nM,1);
plm.Dset = cell(plm.nM,1);
plm.nC   = zeros(plm.nM,1);
plm.nD   = zeros(plm.nM,1);
if Nt || Nf,
    
    % Load FSL style t contrasts
    tcon = cell(Nt,1);
    for t = 1:Nt,
        tmp = palm_miscread(opts.t{t},[],[],opts.precision);
        if any(strcmp(tmp.readwith,{'vestread','csvread','load'})),
            tcon{t} = tmp.data;
        else
            error('Invalid t contrast file: %s',opts.t{t});
        end
    end
    
    % Load FSL style F contrasts
    fcon = cell(Nt,1);
    for t = 1:Nt,
        if ~ isempty(opts.f{t}),
            tmp = palm_miscread(opts.f{t},[],[],opts.precision);
            if any(strcmp(tmp.readwith,{'vestread','csvread','load'})),
                fcon{t} = tmp.data;
            else
                error('Invalid F contrast file: %s',opts.f{t});
            end
        end
    end
    
    % For each valid design, assemble the contrasts.
    for m = 1:plm.nM,
        if Nt == 1;
            t = 1;
        else
            t = m;
        end
        c = 1;
        for j = 1:size(tcon{t},1),
            plm.Cset{m}{c} = tcon{t}(j,:)';
            c = c + 1;
        end
        for j = 1:size(fcon{t},1),
            if ~ isempty(fcon{t}),
                plm.Cset{m}{c} = tcon{t}(logical(fcon{t}(j,:)),:)';
                c = c + 1;
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        for c = 1:plm.nC(m),
            plm.Dset{m}{c} = eye(plm.nY);
        end
        plm.nD(m) = numel(plm.Dset{m});
    end
    
elseif Ncon,
    
    % Load MSET style contrasts (all contrast pairs)
    Ccon = cell(Ncon,1);
    Dcon = cell(Ncon,1);
    for con = 1:Ncon,
        tmp = palm_miscread(opts.Ccon{con},[],[],opts.precision);
        if strcmpi(tmp.readwith,'mset'),
            Ccon{con} = tmp.data;
        else
            error(['Files given to the option "-con" must be in .mset format.\n' ...
                'For .csv/.con/.fts files, use "-t" or "-f".%s'],'');
        end
        if isempty(opts.Dcon{con}),
            for c = 1:numel(Ccon{con}),
                Dcon{con}{c} = eye(plm.nY);
            end
        else
            tmp = palm_miscread(opts.Dcon{con},[],[],opts.precision);
            if strcmpi(tmp.readwith,'mset'),
                Dcon{con} = tmp.data;
            else
                error(['Files given to the option "-con" must be in .mset format.\n' ...
                    'For .csv/.con/.fts files, use "-t" or "-f".%s'],'');
            end
        end
    end
    
    % Assign the contrast sets to the design matrices
    for m = 1:plm.nM,
        if Ncon == 1;
            con = 1;
        else
            con = m;
        end
        plm.Cset{m} = Ccon{con};
        plm.Dset{m} = Dcon{con};
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
        for c = 1:plm.nC(m),
            plm.Cset{m}{c} = plm.Cset{m}{c}';
        end
        for d = 1:plm.nD(m),
            plm.Dset{m}{d} = plm.Dset{m}{d}';
        end
    end
else
    % If no constrasts were at all specified:
    for m = 1:plm.nM,
        if size(plm.Mset{m},2) == 1,
            
            % If there is only 1 regressor, test its effect both
            % positive and negative.
            % The statistic will be t or v, depending on the number of VGs.
            plm.Cset{m}{1} = 1;
            plm.Cset{m}{2} = -1;
            plm.Dset{m}{1} = eye(plm.nY);
            plm.Dset{m}{2} = eye(plm.nY);
        else
            % Otherwise, run an F-test over all regressors in the design matrix.
            % The statistic will be F or G, depending on the number of VGs.
            plm.Cset{m}{1} = eye(size(plm.Mset{m},2));
            plm.Dset{m}{1} = eye(plm.nY);
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
end

% Ranks of the contrasts
plm.rC = cell(plm.nM,1);
plm.rD = plm.rC;
for m = 1:plm.nM,
    plm.rC{m} = zeros(plm.nC(m),1);
    plm.rD{m} = plm.rC{m};
    for c = 1:plm.nC(m),
        plm.rC{m}(c) = rank(plm.Cset{m}{c});
        plm.rD{m}(c) = rank(plm.Dset{m}{c});
    end
end
plm.rC0 = plm.rC; % the rC can be changed for z-scores, but not rC0.

% If only the t or F tests are to be performed
if opts.tonly,
    for m = 1:plm.nM,
        for c = plm.nC(m):-1:1,
            if plm.rC{m}(c) > 1,
                plm.Cset{m}(c) = [];
                plm.Dset{m}(c) = [];
                plm.rC{m}(c)   = [];
                plm.rC0{m}(c)   = [];
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
elseif opts.fonly,
    for m = 1:plm.nM,
        for c = plm.nC(m):-1:1,
            if plm.rC{m}(c) == 1,
                plm.Cset{m}(c) = [];
                plm.Dset{m}(c) = [];
                plm.rC{m}(c)   = [];
                plm.rC0{m}(c)   = [];
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
end

% Some more sanity checks
for m = 1:plm.nM,
    for c = 1:plm.nC(m),
        if any(isnan(plm.Cset{m}{c}(:))) || any(isinf(plm.Cset{m}{c}(:))),
            error('The constrasts cannot contain NaN or Inf.');
        end
        if size(plm.Cset{m}{c},1) ~= size(plm.Mset{m},2),
            error('The size of one or more contrasts don''t match the size of the respective design matrix.')
        end
    end
    for c = 1:plm.nD(m),
        if any(isnan(plm.Dset{m}{c}(:))) || any(isinf(plm.Dset{m}{c}(:))),
            error('The constrasts cannot contain NaN or Inf.');
        end
    end
end
if opts.MV
    if any(plm.nC ~= plm.nD),
        error('The number of C and D contrasts must be the same');
    end
    for m = 1:plm.nM,
        for d = 1:plm.nD(m),
            if size(plm.Dset{m}{d},1) ~= plm.nY,
                error('The size of one or more MV contrasts don''t match the number of modalities.');
            end
        end
    end
end
for m = 1:plm.nM,
    for c = 1:plm.nC(m),
        if opts.concordant && plm.rC{m}(c) > 1,
            error(['Cannot use the "-concordant" option with F-tests (inadmissible).\n'...
                'Use "-tonly" to run just the t-tests, or remove "-concordant".%s'],'');
        end
    end
end

% Check if the contrasts have all the same rank for correction over
% contrasts. If not, convert to zstat.
if opts.corrcon && ~ opts.zstat,
    rC1 = plm.rC{1}(1);
    rD1 = plm.rD{1}(1);
    for m = 1:plm.nM,
        brflag = false;
        for c = 1:plm.nC(m),
            if rC1 ~= plm.rC{m}(c) || rD1 ~= plm.rD{m}(c);
                warning([...
                    'Not all contrasts have the same rank, and the option "-corrcon" was used.\n' ...
                    '         Adding the option "-zstat" automatically.%s'],'');
                opts.zstat = true;
                brflag = true;
                break;
            end
        end
        if brflag,
            break;
        end
    end
end

% Partition the model according to the contrasts and design matrix.
% The partitioning needs to be done now, because of the need for
% synchronised permutations/sign-flips
if ~ opts.cmcx,
    seqtmp = zeros(plm.N,sum(plm.nC));
    j = 1;
    plm.seq = cell(plm.nM,1);
    for m = 1:plm.nM,
        plm.seq{m} = cell(plm.nC(m),1);
        for c = 1:plm.nC(m),
            Xtmp = palm_partition(plm.Mset{m}(:,:,1),plm.Cset{m}{c},opts.pmethodp);
            [~,~,plm.seq{m}{c}] = unique(Xtmp,'rows');
            seqtmp(:,j) = plm.seq{m}{c};
            j = j + 1;
        end
    end
    tmp = sum(diff(seqtmp,1,2).^2,2);
    if (opts.corrcon || opts.npccon || opts.syncperms) && any(tmp(:) ~= 0),
        warning([ ...
            'You chose to correct over contrasts, or run NPC\n'             ...
            '         between contrasts, but with the design(s) and,\n'     ...
            '         contrasts given it is not possible to run\n'          ...
            '         synchronised permutations without ignoring repeated\n'...
            '         elements in the design matrix (or matrices). To\n'    ...
            '         solve this, adding the option "-cmcx" automatically.%s\n'],'');
        opts.cmcx = true;
    end
end
if opts.cmcx,
    plm.seq = cell(plm.nM,1);
    for m = 1:plm.nM,
        plm.seq{m} = cell(plm.nC(m),1);
        for c = 1:plm.nC(m),
            plm.seq{m}{c} = (1:plm.N)';
        end
    end
    if opts.corrcon || opts.npccon,
        opts.syncperms = true;
    end
end

% Make sure not too many components are asked if CCA or PLS are used
if opts.CCA || opts.PLS,
    if opts.ccaorplsparm > plm.nY,
        error(['Cannot ask more canonical correlations (CCA) or \n', ...
            'score vectors (PLS) (k=%d) than the number of modalities (#(-i)=%d).\n'],...
            opts.ccaorplsparm,plm.nY);
    end
    for m = 1:plm.nM,
        for c = 1:plm.nC(m),
            if opts.ccaorplsparm > plm.rC{m}(c),
                error(['Cannot ask more canonical correlations (for CCA) or score \n', ...
                    'vectors (for PLS) than the rank of the contrast (k=%d > rank=%d).\n', ...
                    'Please check design %d, contrast %d.'],opts.ccaorplsparm,plm.rC{m}(c),m,c);
            end
        end
    end
end

% Read the exchangeability blocks. If none is specified, all observations
% are assumed to be in the same large block. Also treat the legacy format of
% a single column for the EBs.
if isempty(opts.eb),
    plm.EB = [];
    if opts.within || opts.whole,
        error([ ...
            'Options -within and/or -whole require a file defining\n' ...
            '         the exchangeability blocks (option -eb).\n%s'],'');
    end
else
    plm.EB = palm_miscread(opts.eb);
    plm.EB = plm.EB.data;
    if isvector(plm.EB),
        if opts.within && opts.whole, % within + whole block shuffling
            plm.EB = [+ones(plm.N,1) +plm.EB(:) (1:plm.N)'];
        elseif opts.whole             % whole-block shuffling
            plm.EB = [+ones(plm.N,1) -plm.EB(:) (1:plm.N)'];
        else                          % within-block shuffling (this is the default, and not meant to be changed)
            plm.EB = [-ones(plm.N,1) +plm.EB(:) (1:plm.N)'];
        end
    elseif opts.within || opts.whole,
        warning([ ...
            'Options -within and/or -whole ignored, as the file defining\n' ...
            '         the exchangeability blocks (option -eb) already \n' ...
            '         defines how the data should be shuffled.%s'],'');
    elseif ~ isempty(plm.subjidx),
        error('Cannot use "-subjidx" with multi-level blocks.');
    end
    plm.EB = palm_reindex(plm.EB,'fixleaves');
end

% Load/define the variance groups.
if opts.singlevg,
    % If single VG, it's all ones
    plm.VG = ones(plm.N,1);
elseif strcmpi(opts.vg,'auto'),
    if isempty(plm.EB),
        % If auto, but there are no exchangeability blocks, it's all ones too
        plm.VG = ones(plm.N,1);
    else
        % Generate an initial dependence tree, to be used to define variance groups.
        % The tree used for the permutations later require the design matrix, and
        % varies for each contrast -- all to be taken care of later.
        Ptree  = palm_tree(plm.EB,(1:plm.N)');
        plm.VG = palm_ptree2vg(Ptree);
    end
else
    % The automatic variance groups can be overriden if the user specified
    % a file with the custom definitions.
    plm.VG = palm_miscread(opts.vg);
    plm.VG = plm.VG.data;
end
if ~ isempty(plm.subjidx) && size(plm.VG,1) ~= plm.N,
    plm.VG = plm.VG(plm.subjidx,:);
end
[tmp,~,plm.VG] = unique(plm.VG);
plm.nVG = numel(tmp);
if plm.nVG == 1, opts.singlevg = true; end

% MV can't yet be used if nVG>1, although NPC remains an option
if opts.MV && plm.nVG > 1,
    error('There are more than one variance group. MV cannot be used (but NPC can).');
end

% There should be no more missing indicators than modalities, or designs.
% These need to be either 1 or the same as the corresponding numbers of
% modalities/designs.
if Nimiss > Ni,
    error([...
        'There are more missing indicators supplied with "-imiss" (%d) than\n'...
        'modalities supplied with "-i" (%d)'],Nimiss,Ni);
elseif Nimiss > 1 && Nimiss ~= Ni,
    error([...
        'The number of missing indicators supplied with "-imiss" (%d) is larger,\n'...
        'than 1, but still not the same as the number of modalities supplied with\n'...
        'the option "-i" (%d).'],Nimiss,Ni);
end
if Ndmiss > Nd,
    error([...
        'There are more missing indicators supplied with "-dmiss" (%d) than\n'...
        'designs supplied with "-d" (%d)'],Nimiss,Ni);
elseif Ndmiss > 1 && Ndmiss ~= Nd,
    error([...
        'The number of missing indicators supplied with "-dmiss" (%d) is larger,\n'...
        'than 1, but still not the same as the number of modalities supplied with\n'...
        'the option "-d" (%d).'],Ndmiss,Nd);
end

% Load the missing indicators for the data ("imiss"):
for i = 1:Nimiss,
    if strcmpi(opts.imiss{i},'none'),
        if isempty(plm.subjidx)
            tmp = zeros(plm.N,1);
        else
            tmp = zeros(size(plm.subjidx,1),1);
        end
    else
        tmp = palm_miscread(opts.imiss{i});
        tmpfname = tmp.filename;
        tmp = tmp.data;
        if ~ isempty(plm.subjidx) && size(tmp,1) ~= plm.N,
            tmp = tmp(plm.subjidx,:);
        end
        checkmiss(tmp,tmpfname,plm.N);
    end
    plm.Ymiss{i} = tmp;
end
if Nimiss == 1,
    for i = 2:Ni,
        plm.Ymiss{i} = plm.Ymiss{1};
    end
end

% Load the missing indicators for the design ("dmiss"):
for d = 1:Ndmiss,
    if strcmpi(opts.dmiss{d},'none'),
        if isempty(plm.subjidx)
            tmp = zeros(plm.N,1);
        else
            tmp = zeros(size(plm.subjidx,1),1);
        end
    else
        tmp = palm_miscread(opts.dmiss{d});
        tmpfname = tmp.filename;
        tmp = tmp.data;
        if ~ isempty(plm.subjidx) && size(tmp,1) ~= plm.N,
            tmp = tmp(plm.subjidx,:);
        end
        checkmiss(tmp,tmpfname,plm.N);
    end
    plm.Mmiss{d} = tmp;
end
if Ndmiss == 1,
    for d = 2:Nd,
        plm.Mmiss{d} = plm.Mmiss{1};
    end
end
for d = 1:Ndmiss,
    if any(size(plm.Mmiss{d}) ~= size(plm.Mset{d})),
        if strcmpi(opts.dmiss{d},'none'),
            plm.Mmiss{d} = repmat(plm.Mmiss{d},[1 size(plm.Mset{d},2)]);
        else
            error([ ...
                'The missing data indicator ("-dmiss") must have\n', ...
                'the same size as the respective design.%s'],'');
        end
    end
end

% If only data or design missing indicators are missing, fill the other
% with all-false indicators.
if     Nimiss && ~ Ndmiss,
    for m = 1:plm.nM,
        plm.Mmiss{m} = false(size(plm.Mset{m}));
    end
elseif Ndmiss && ~ Nimiss
    for y = 1:plm.nY,
        plm.Ymiss{y} = false(size(plm.Yset{y}));
    end
end

% Remove the variance groups with tiny sample sizes?
if plm.nVG > 1 && ~ opts.removevgbysize && (opts.vgdemean || opts.ev4vg) && ...
        any(sum(bsxfun(@eq,plm.VG,unique(plm.VG)'),1) == 1),
    warning([...
        'The options "-vgdemean" and "-ev4vg" require that observations\n' ...
        '         in variance groups of size 1 are removed.\n' ...
        '         Enabling the option "-removevgbysize 1"%s.'],'');
    opts.removevgbysize = 1;
end
if ~ opts.removevgbysize,
    tmpwarned = false;
    for u = 1:plm.nVG,
        if sum((plm.VG == u),1) == 1,
            if ~ tmpwarned,
                warning([...
                    'There are variance groups with just one observation.\n' ...
                    '         Consider using the option "-removevgbysize 1" to improve the\n' ...
                    '         variance estimates (at the cost of reducing sample size).%s'],'');
                tmpwarned = true;
            end
        end
    end
end
if opts.removevgbysize > 0,
    
    % Indices of the observations to keep
    uVG = unique(plm.VG)';
    idxvg = sum(bsxfun(@eq,plm.VG,uVG),1) <= opts.removevgbysize;
    idx   = any(bsxfun(@eq,plm.VG,uVG(~idxvg)),2);
    
    % Modify all data as needed
    for y = 1:plm.nY,
        plm.Yset{y} = plm.Yset{y}(idx,:);
    end
    if ~ isempty(plm.EB),
        plm.EB = plm.EB(idx,:);
    end
    for m = 1:plm.nM,
        plm.Mset{m} = plm.Mset{m}(idx,:);
    end
    plm.N = sum(idx);
    [tmp,~,plm.VG] = unique(plm.VG(idx));
    plm.nVG = numel(tmp);
end

% Add one regressor for each variance group if requested
if opts.ev4vg,
    for m = 1:plm.nM,
        Mvg = zeros(plm.N,plm.nVG);
        V = unique(plm.VG);
        for v = 1:plm.nVG,
            Mvg(plm.VG == V(v),v) = 1;
        end
        rM   = round(sum(diag(plm.Mset{m}*pinv(plm.Mset{m}))));
        Mnew = horzcat(plm.Mset{m},Mvg);
        if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG),
            plm.Mset{m} = Mnew;
            nadded      = plm.nVG;
        else
            Mnew = Mnew(:,1:end-1);
            if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG - 1),
                plm.Mset{m} = Mnew;
                nadded      = plm.nVG - 1;
            else
                error([ ...
                    'It was not possible to add one regressor for each variance group\n' ...
                    'perhaps because they already exist in the design. Check your design\n' ...
                    'matrix and maybe consider including these regressors manually.%s'],'');
            end
        end
        for c = 1:plm.nC(m),
            plm.Cset{m}{c} = vertcat(plm.Cset{m}{c},...
                zeros(nadded,size(plm.Cset{m}{c},2)));
        end
    end
end

% Remove intercept from the design for the options -demean and -vgdemean
if opts.demean || opts.vgdemean,
    for m = 1:plm.nM,
        tmp = plm.Mset{m};
        siz = size(tmp);
        intercp = all(bsxfun(@eq,reshape(plm.Mset{m}(1,:),[1 siz(2:end)]),plm.Mset{m}),1);
        if any(intercp),
            for c = 1:plm.nC(m),
                if any(intercp*plm.Cset{m}{c}~=0,2),
                    error([ ...
                        'Contrast %d (and perhaps others) tests the intercept. This means\n' ...
                        'that the options "-demean" and "-vgdemean" cannot be used.\n' ...
                        'If "-demean" was added to calculate Pearson''s "r" or the "R^2"\n' ...
                        'note that these statistics cannot be computed for constant variables.%s'],c,'');
                else
                    plm.Cset{m}{c}(intercp,:) = [];
                end
            end
            plm.Mset{m}(:,intercp) = [];
        end
    end
end

% Mean center data and design (-demean)
if opts.demean,
    for m = 1:plm.nM,
        plm.Mset{m} = bsxfun(@minus,plm.Mset{m},mean(plm.Mset{m},1));
    end
    for y = 1:plm.nY,
        plm.Yset{y} = bsxfun(@minus,plm.Yset{y},mean(plm.Yset{y},1));
    end
end

% Mean center data and design, within VG
if opts.vgdemean,
    
    % For each VG
    V = unique(plm.VG);
    for v = 1:plm.nVG,
        vidx = plm.VG == V(v);
        
        % Demean design within VG
        for m = 1:plm.nM,
            plm.Mset{m}(vidx,:) = bsxfun(@minus,...
                plm.Mset{m}(vidx,:),mean(plm.Mset{m}(vidx,:),1));
        end
        
        % Demean data within VG
        for y = 1:plm.nY,
            plm.Yset{y}(vidx,:) = bsxfun(@minus,...
                plm.Yset{y}(vidx,:),mean(plm.Yset{y}(vidx,:),1));
        end
    end
end

% Number of tests to be selected for the low rank approximation
if opts.accel.lowrank,
    if plm.nVG > 1,
        error('The option "-accel lowrank" cannot be used with more than one variance group.');
    end
    if opts.nP0 == 0,
        error('With lowrank approximation you must indicate a larger-than-zero number of permutations.');
    end
    if opts.nP0 < plm.N*(plm.N+1)/2,
        error([ ...
            'Too few permutations selected to use with lowrank approximation.\n' ...
            'Use at least N*(N+1)/2 = %d to note a speed difference and have reasonably accurate results.\n'...
            'Otherwise, don''t bother using lowrank approximation.\n'],plm.N*(plm.N+1)/2);
    end
    if opts.spatial.do,
        warning([ ...
            'There isn''t much benefit in using lowrank approximation with spatial statistics\n' ...
            '         like TFCE and cluster extent and/or mass. These cannot be accelerated with this\n' ...
            '         method, and the overall gain will be minimal. Consider other approximation methods,\n' ...
            '         or run the full permutation test, or just drop spatial statistics.%s'],'');
    end
    plm.nsel = zeros(plm.nY,1);
    if opts.accel.lowrank_val <= 1,
        for y = 1:plm.nY,
            plm.nsel(y) = ceil(opts.accel.lowrank_val*plm.Ysiz(y));
        end
    elseif opts.accel.lowrank_val > 1
        plm.nsel(1:end) = ceil(opts.accel.lowrank_val);
    elseif isnan(opts.accel.lowrank_val),
        plm.nsel(1:end) = plm.N*(plm.N+1)/2;
    end
end

% ==============================================================
function checkmiss(A,Afname,N)
% Check if the missing data indicators are sane.
for a = 1:size(A,2),
    U = unique(A(:,a));
    if size(A,1) ~= N,
        error([ ...
            'The missing data indicators ("-imiss" and "-dmiss") must have\n', ...
            'the same number of observations as the data and design.'],'');
    elseif ...
            (numel(U) >  2) || ...
            (numel(U) == 2 && ~ any(U == 0)) || ...
            (numel(U) == 2 && ~ any(U(U~=0) == [-1 1 2])) || ...
            (numel(U) == 1 && U ~= 0),
        error([ ...
            'The missing data indicators ("-imiss" and "-dmiss") must have\n', ...
            'no more than two unique values per column, one being 0, the\n', ...
            'the other being either -1, 1, or 2.\n', ...
            'Consult the documentation for details.\n'...
            '- Filename: %s\n',...
            '- Column: %d (possibly also others)'],Afname,a);
    end
end
