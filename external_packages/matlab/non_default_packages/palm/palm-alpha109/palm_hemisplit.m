function palm_hemisplit(varargin)
% Split back bh into lh and rh for files that have previously
% been merged with palm_hemimerge.
% Output file names have prefix "lh" and "rh".
%
% palm_hemisplit [-v #vtx] [-f #fac] <files>
%
% Wildcards are accepted.
%
% -v <nV> : Specify the number of vertices for the left hemisphere.
% -f <nF> : Specify the number of faces for the left hemisphere.
% 
% For surfaces (meshes) and curvature files, both -v and -f are necessary.
% For dpv or dpx, use -v for the same purpose. For dpf, use -f. 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jun/2016
% http://brainder.org

% List of files (with wildcards)
a   = 1;
j   = 1;
nVL = [];
nFL = [];
while a <= nargin,
    if a < nargin && strcmpi(varargin{a},'-v'),
        nVL = varargin{a+1};
        if ischar(nVL),
            nVL = str2double(nVL);
        end
        a = a + 2;
    elseif a < nargin && strcmpi(varargin{a},'-f'),
        nFL = varargin{a+1};
        if ischar(nFL),
            nFL = str2double(nFL);
        end
        a = a + 2;
    else
        F = dir(varargin{a});
        for f = numel(F):-1:1,
            if F(f).name(1) == '.',
                F(f) = [];
            else
                Flist{j} = F(f).name;
                j = j + 1;
            end
        end
        a = a + 1;
    end
end
Flist = flipud(Flist');

% For each input file.
for f = 1:numel(Flist),
    fprintf('Working on: %s\n',Flist{f});
    
    B = palm_miscread(Flist{f});
    L = B; R = B;
    switch B.readwith,
        case {'load','csvread','fs_load_mgh'},
            
            nVB = size(B.data,1);
            if isempty(nVL),
                nVL = nVB/2;
            end
            L.data = B.data(1:nVL,:,:,:);
            R.data = B.data(nVL+1:end,:,:,:);
            
        case 'dpxread',
            
            nXB = size(B.data,1);
            if strcmpi(Flist{f}(end-2:end),'dpf'),
                if isempty(nFL),
                    nXL = nXB/2;
                else
                    nXL = nFL;
                end
            else
                if isempty(nVL),
                    nXL = nXB/2;
                else
                    nXL = nVL;
                end
            end
            L.data      = B.data(1:nXL,:,:,:);
            R.data      = B.data(nXL+1:end,:,:,:);
            L.extra.crd = B.extra.crd(1:nXL,:,:,:);
            R.extra.crd = B.extra.crd(nXL+1:end,:,:,:);
            L.extra.idx = (1:size(L.data,1))';
            R.extra.idx = (1:size(R.data,1))';
            
        case {'srfread','fs_read_surf'},
            
            nVB = size(B.data.vtx,1);
            nFB = size(B.data.fac,1);
            if xor(isempty(nVL),isempty(nFL)),
                warning(...
                    ['With surfaces (meshes), either both "-v" and "-f" must be supplied, or neither.\n'...
                    'Skipping: %s.'],Flist{f});
            else
                if isempty(nVL) && isempty(nFL),
                    nVL = nVB/2;
                    nFL = nFB/2;
                end
                L.data.vtx = B.data.vtx(1:nVL,:);
                R.data.vtx = B.data.vtx(nVL+1:end,:);
                L.data.fac = B.data.fac(1:nFL,:);
                R.data.fac = B.data.fac(nFL+1:end,:)-nVL;
            end
            
        case 'fs_read_curv',
            
            nVB = size(B.data,1);
            nFB = B.extra.fnum;
            if xor(isempty(nVL),isempty(nFL)),
                warning(...
                    ['With curvatures, either both "-v" and "-f" must be supplied, or neither.\n'...
                    'Skipping: %s.'],Flist{f});
            else
                if isempty(nVL) && isempty(nFL),
                    nVL = nVB/2;
                    nFL = nFB/2;
                end
                L.data      = B.data(1:nVL,:,:,:);
                R.data      = B.data(nVL+1:end,:,:,:);
                L.extra.fnum = nFL;
                R.extra.fnum = nFB-nFL;
            end
            
        otherwise
            warning('Cannot deal with files read with %s. Skipping: %s\n',B.readwith,Flist{f});
    end
    
    % Adjust filenames and save
    filename = Flist{f};
    if any(strcmpi(filename(end-3:end),{'.mgh','.mgz'})),
        filename = filename(1:end-4);
    end
    if strcmpi(Flist{f}(1:2),'bh'),
        L.filename = filename;
        L.filename(1:2) = 'lh';
        R.filename = filename;
        R.filename(1:2) = 'rh';
    else
        L.filename = strcat('lh_',filename);
        R.filename = strcat('rh_',filename);
    end
    palm_miscwrite(L);
    palm_miscwrite(R);
end
