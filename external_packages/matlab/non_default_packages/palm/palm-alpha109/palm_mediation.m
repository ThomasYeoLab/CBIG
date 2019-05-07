function palm_mediation(varargin)
% Perform a mediation analysis using the strategy proposed by
% Baron & Kenny (1986).
%
% palm_mediation -indep <file_independent_variable> ... 
%                -dep   <file_dependent_variable>   ...
%                -med   <file_mediator_variable>    ...
%                -o     <output_prefix>             ...
%                [ -m   <file_mask> ] [other options]
% 
% Type "palm_mediation" without arguments for a description
% of these options.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2015
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

% Remove trailing empty arguments. This is useful for some Octave versions.
while numel(varargin) > 0 && isempty(varargin{1}),
    varargin(1) = [];
end
nargin = numel(varargin);
if nargin == 0,
    palm_help('logo');
    fprintf('                        * Mediation Analysis *\n');
    fprintf('=======================================================================\n\n');
    fprintf('The main options are:\n\n');
    fprintf('-indep <file> : File with the independent variable.\n\n');
    fprintf('-dep <file> : File with the dependent variable.\n\n');
    fprintf('-med <file> : File with the putative mediator.\n\n');
    fprintf('-m <file> : Mask.\n\n');
    fprintf('-o <prefix> : Output prefix.\n\n');
    palm_help('date');
    return;
end

% Load the defaults
medi = palm_defaults;
medi.transposedata = false;
medi.fdr = false;

% List of forbidden options with mediation analysis:
forbidden_opts = {                                    ...
    '-i', '-transposedata',                           ...
    '-t', '-f', '-fonly', '-con', '-conskipcount',    ...
    '-save1-p', '-nounivariate', '-verbosefilenames', ...
    '-npc', '-npcmod', '-npccon',                     ...
    '-mv', '-inputmv', '-pearson', '-noranktest',     ...
    '-C', '-Cuni', '-Cmv', '-Cnpc', '-Cstat',         ...
    '-corrmod', '-corrcon'};

% List of mandatory options (even if the user doesn't specify,
% these will be there:
palmopts = {           ...
    '-designperinput', ...
    '-twotail',        ...
    '-logp',           ...
    '-saveglm',        ...
    '-savemask',       ...
    '-forceintersectionmask'};

% Take the input arguments
medi.nuisance = [];
medi.i = cell(3,1);
Nevd = sum(strcmp(varargin,'-evperdat')); % number of EV per datum inputs
a = 1; ev = 1;
while a <= nargin,
    switch varargin{a},
        case {'-indep','-independent'},
            
            % Get the filename for the independent variable (X).
            medi.i{1} = varargin{a+1};
            a = a + 2;
            
        case {'-med','-mediator'},
            
            % Get the filename for the mediator (M).
            medi.i{2} = varargin{a+1};
            a = a + 2;
            
        case {'-dep','-dependent'},
            
            % Get the filenames for the dependent variable (Y).
            medi.i{3} = varargin{a+1};
            a = a + 2;
            
        case {'-nui','-nuisance'},
            
            % Get the filename for the nuisance variables (Z).
            medi.nuisance = varargin{a+1};
            a = a + 2;
            
        case '-evperdat',
            
            % Use one EV per datum?
            medi.evperdat = true;
            medi.evdatfile{ev} = varargin{a+1};
            ev = ev + 1;
            a  = a + 2;
            
        case '-o',
            
            % Output prefix for the files to be saved.
            medi.o = varargin{a+1};
            a = a + 2;
            
        case '-noniiclass',
            medi.useniiclass = false;
            palmopts{end+1} = varargin{a};
            a = a + 1;
            
        case '-fdr',
            medi.fdr = true;
            palmopts{end+1} = varargin{a};
            a = a + 1;
            
        case forbidden_opts,
            errstr = sprintf('\t%s\n',forbidden_opts{:});
            error('The options below cannot be used with palm_mediation:\n%s',errstr);
            
        otherwise
            palmopts{end+1} = varargin{a};
            a = a + 1;
    end
end

% ==============================================================
% Stage I: Prepare the 3 GLMs.
% ==============================================================
fprintf('=======================================================================\n');
fprintf('Mediation stage I: Data and model preparation.\n');

% Read input files:
N = zeros(3,1);
for i = 1:3,
    fprintf('Reading input %d/%d: %s\n',i,3,medi.i{i});
    [I{i},maskstruct{i},~,~,~,Itmp{i}] = palm_ready(medi.i{i},[],medi);
    [siz(i,1),siz(i,2),siz(i,3)] = size(maskstruct{i}.data);
    N(i) = size(I{i},1);
end

% Check sizes:
tmp = size(unique(siz,'rows'),1);
masksiz = prod(siz,2);
if ~ (numel(unique(N)) == 1 && ...
        (tmp == 1 || ...
        (tmp == 2 && min(masksiz) == 1))),
    error('Input data of incompatible sizes.');
end
N = N(1);

% Intersection mask:
maskinter = [];
for i = 1:3,
    if numel(maskstruct{i}.data) > 1,
        r = i;
        if isempty(maskinter),
            maskinter = maskstruct{i}.data;
        else
            maskinter = maskinter & maskstruct{i}.data;
        end
    end
end

% If none of the main inputs is voxelwise, but there are nuisance voxelwise
% EVs, make the mask from the first of these:
if isempty(maskinter) && Nevd > 0,
    r = 0;
    fprintf('Reading file: %s\n',medi.evdatfile{1});
    [~,maskstruct,~,~,~,EVtmp] = palm_ready(medi.evdatfile{1},[],medi);
    maskinter = maskstruct.data;
end

% Expansion of M: If independent or nuisance is voxelwise, but not the
% mediator, expand it to be voxelwise too:
if (prod(siz(1,:)) > 1 && prod(siz(2,:)) ==  1) || Nevd > 1,
    fprintf('Expanding putative mediator to match the size of inputs.\n');
    if r == 0,
        Inew = EVtmp;
    else
        Inew = Itmp{r};
    end
    Inew.data = bsxfun(@times,maskinter,permute(I{2},[2 3 4 1]));
    
    % File name of the temporary file:
    [~,ifnam,ifext] = fileparts(Itmp{2}.filename);
    if strcmpi(ifext,'.gz'),
        [~,ifnam,~] = fileparts(ifnam);
    end
    [~,rfnam,rfext] = fileparts(Inew.filename);
    if strcmpi(rfext,'.gz'),
        [~,~,rfext2] = fileparts(rfnam);
        rfext = strcat(rfext2,rfext);
    end
    Inew.filename = strcat(medi.o,'_expanded_',ifnam,rfext);
    
    % Save it
    palm_miscwrite(Inew);
    medi.i{2} = Inew.filename;
    Mexpanded = true;
else
    Mexpanded = false;
end

% Expansion of Y: If independent or nuisance is voxelwise, but not the
% dependent, expand it to be voxelwise too:
if (prod(siz(1,:)) > 1 && prod(siz(3,:)) ==  1) || Nevd > 1,
    fprintf('Expanding dependent variable to match the size of inputs.\n');
    if r == 0,
        Inew = EVtmp;
    else
        Inew = Itmp{r};
    end
    Inew.data = bsxfun(@times,maskinter,permute(I{3},[2 3 4 1]));
    
    % File name of the temporary file:
    [~,ifnam,ifext] = fileparts(Itmp{3}.filename);
    if strcmpi(ifext,'.gz'),
        [~,ifnam,~] = fileparts(ifnam);
    end
    [~,rfnam,rfext] = fileparts(Inew.filename);
    if strcmpi(rfext,'.gz'),
        [~,~,rfext2] = fileparts(rfnam);
        rfext = strcat(rfext2,rfext);
    end
    Inew.filename = strcat(medi.o,'_expanded_',ifnam,rfext);
    
    % Save it
    palm_miscwrite(Inew);
    medi.i{3} = Inew.filename;
    Yexpanded = true;
else
    Yexpanded = false;
end
clear('I','maskstruct','Itmp','EVtmp','Inew')

% Load the nuisance variables:
if isempty(medi.nuisance),
    Z = ones(N,1);
else
    Z = palm_miscread(medi.nuisance);
    Z = Z.data;
    if ndims(Z) ~= 2,
        error([...
            'THe variable with nuisance must be a vector or matrix.\n',...
            'For voxelwise nuisance, use "-evperdat"\n.%s'],'');
    end
end

% Add dummy EVs for the nuisance -evperdat:
if Nevd > 0,
    Z = horzcat(zeros(N,Nevd),Z);
end

% Prepare options for the GLM #1:
fprintf('Preparing options of PALM.\n');
opts_step1 = { ...
    '-i',medi.i{3}, ...
    '-d',fullfile(sprintf('%s_glm1+2.mat',medi.o)), ...
    '-t',fullfile(sprintf('%s_glm1+2.con',medi.o)) };
if masksiz(1) == 1,
    X = palm_miscread(medi.i{1},medi.useniiclass,medi.o);
    X = X.data;
elseif masksiz(1) > 1,
    X = zeros(N,1);
    opts_step1(end+1:end+4) = { ...
        '-evperdat',medi.i{1},size(Z,2)+1,1};
end
palm_vestwrite(fullfile(sprintf('%s_glm1+2.mat',medi.o)),[Z X]);
palm_vestwrite(fullfile(sprintf('%s_glm1+2.con',medi.o)),[zeros(1,size(Z,2)) 1]);
for ev = 1:Nevd,
    opts_step1(end+1:end+4) = { ...
        '-evperdat',medi.evdatfile{ev},ev,1};
end

% Prepare options for the GLM #2:
opts_step2 = opts_step1;
opts_step2{find(strcmp('-i',opts_step2))+1} = medi.i{2};
for f = find(strcmp('-evperdat',opts_step2)),
    opts_step2{f+3} = 2;
end

% Prepare options for the GLM #3a:
opts_step3a = { ...
    '-i',medi.i{3}, ...
    '-d',fullfile(sprintf('%s_glm3.mat',medi.o)), ...
    '-t',fullfile(sprintf('%s_glm3.con',medi.o)) };
if masksiz(1) == 1,
    X = palm_miscread(medi.i{1},medi.useniiclass,medi.o);
    X = X.data;
elseif masksiz(1) > 1,
    X = zeros(N,1);
    opts_step3a(end+1:end+4) = { ...
        '-evperdat',medi.i{1},size(Z,2)+1,3};
end
if masksiz(2) == 1 && ~ Mexpanded,
    M = palm_miscread(medi.i{2},medi.useniiclass,medi.o);
    M = M.data;
else
    M = zeros(N,1);
    opts_step3a(end+1:end+4) = { ...
        '-evperdat',medi.i{2},size(Z,2)+2,3};
end
palm_vestwrite(fullfile(sprintf('%s_glm3.mat',medi.o)),[Z X M]);
palm_vestwrite(fullfile(sprintf('%s_glm3.con',medi.o)),[zeros(1,size(Z,2)+1) 1]);
for ev = 1:Nevd,
    opts_step3a(end+1:end+4) = { ...
        '-evperdat',medi.evdatfile{ev},ev,3};
end

% Prepare options for the GLM #3b:
opts_step3b = opts_step3a;
idx = find(strcmpi(opts_step3b,'-evperdat'));
for ev = idx,
    opts_step3b{ev+3} = 1;
end
palmopts3b = palmopts;
idx = find(strcmpi(palmopts3b,'-n'),1,'last');
if isempty(idx),
    palmopts3b(end+1:end+2) = {'-n',1};
else
    palmopts3b{idx+1} = 1;
end
fprintf('Finished stage I.\n');

% ==============================================================
% Stage II: Run PALM.
% ==============================================================
fprintf('=======================================================================\n');
fprintf('Mediation stage II: Permutation test.\n');
palm( ...
    opts_step1{:},opts_step2{:},opts_step3a{:},palmopts{:},...
    '-o',sprintf('%s_1+2+3a',medi.o),'-corrcon','-cmcx');
palm(opts_step3b{:},palmopts3b{:},...
    '-o',sprintf('%s_3b',medi.o));
fprintf('Finished stage II.\n');

% ==============================================================
% Stage III: Check if there is mediation.
% ==============================================================
fprintf('=======================================================================\n');
fprintf('Mediation stage III: Inference.\n');

% Figure out what is the stat name (tstat or vstat) and the file extension:
F = dir(sprintf('%s_1+2+3a_*_fwep_m1_d1.*',medi.o));
kindname = cell(numel(F),1);
statname = cell(numel(F),1);
for f = 1:numel(F),
    F(f).name = F(f).name(numel(medi.o)+9:end);
    [kindname{f},F(f).name] = strtok(F(f).name,'_');
    statname{f} = strtok(F(f).name(2:end),'_');
end
Ykindstr{1} = kindname{find(~strcmpi('tfce',kindname),1)};
if any(strcmpi('tfce',kindname)),
    Ykindstr{2} = 'tfce';
end
statname = statname{1};
[~,fnam,fext] = fileparts(F(1).name);
if strcmpi(fext,'.gz'),
    [~,~,fext2] = fileparts(fnam);
    fext = strcat(fext2,fext);
end

% Check where the COPE in step 1 is larger than in step 3a. This is used
% below to essentially mask the results:
S1 = palm_miscread( ...
    sprintf('%s_1+2+3a_%s_cope_m1_d1%s',medi.o,Ykindstr{1},fext),...
    medi.useniiclass,medi.o);
S3b = palm_miscread( ...
    sprintf('%s_3b_%s_cope%s',medi.o,Ykindstr{1},fext),...
    medi.useniiclass,medi.o);
S3b_mask = double(S1.data) > double(S3b.data);

% Check places where steps 1, 2 and 3a are significant. Do it as a
% conjunction.
for f = 1:numel(Ykindstr),
    S = [];
    for r = 1:3,
        R = palm_miscread(...
            sprintf('%s_1+2+3a_%s_%s_fwep_m%d_d%d%s',medi.o,Ykindstr{f},statname,r,r,fext),...
            medi.useniiclass,medi.o);
        if isempty(S),
            S = double(R.data);
        else
            S = min(S,double(R.data));
        end
    end
    S = S .* S3b_mask;
    S1.data = S;
    S1.filename = sprintf('%s_%s_mediation_fwep%s',medi.o,Ykindstr{f},fext);
    palm_miscwrite(S1);
end

% Repeat the last part if FDR has been choosen:
if medi.fdr,
    for f = 1:numel(Ykindstr),
        S = [];
        for r = 1:3,
            R = palm_miscread(...
                sprintf('%s_1+2+3a_%s_%s_fdrp_m%d_d%d%s',medi.o,Ykindstr{f},statname,r,r,fext),...
                medi.useniiclass,medi.o);
            if isempty(S),
                S = double(R.data);
            else
                S = min(S,double(R.data));
            end
        end
        S = S .* S3b_mask;
        S1.data = S;
        S1.filename = sprintf('%s_%s_mediation_fdrp%s',medi.o,Ykindstr{f},fext);
        palm_miscwrite(S1);
    end
end
fprintf('Finished stage III.\nDone.\n');
