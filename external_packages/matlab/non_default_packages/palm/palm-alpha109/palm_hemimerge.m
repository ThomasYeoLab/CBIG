function palm_hemimerge(varargin)
% Merge data for lh and rh (left and right hemispheres) into a
% single file. Various input files can be given, starting with
% 'lh' and/or 'rh'. The script finds the correct pair, removes
% repeated files, and merges them into files beginning with 'bh'.
%
% palm_hemimerge <files>
%
% Wildcards are accepted.
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jun/2016
% http://brainder.org

% List of files (with wildcards)
j = 1;
for a = 1:nargin,
    F = dir(varargin{a});
    for f = numel(F):-1:1,
        if F(f).name(1) == '.',
            F(f) = [];
        else
            Flist{j} = F(f).name;
            j = j + 1;
        end
    end
end
Flist = flipud(Flist');

% Prepare pairs to merge and output names
F = cell(0,3);
f = 1;
for a = numel(Flist):-1:1,
    if strcmpi(Flist{a}(1:2),'lh'),
        if ~ any(strcmpi(Flist{a},F(:,1))),
            F{f,1} = Flist{a};
            F{f,2} = strcat('rh',Flist{a}(3:end));
            F{f,3} = strcat('bh',Flist{a}(3:end));
            f = f + 1;
        end
    elseif strcmpi(Flist{a}(1:2),'rh'),
        if ~ any(strcmpi(Flist{a},F(:,2))),
            F{f,1} = strcat('lh',Flist{a}(3:end));
            F{f,2} = Flist{a};
            F{f,3} = strcat('bh',Flist{a}(3:end));
            f = f + 1;
        end
    else
        warning('Files for merger must start with "lh" or "rh". Skipping: %s\n',Flist{a});
    end
end

% For each valid pair
for f = 1:size(F,1),
    
    fprintf('Working on: %s and %s\n',F{f,1},F{f,2});
    
    % Load L and R
    L = palm_miscread(F{f,1});
    R = palm_miscread(F{f,2});
    if ~ strcmpi(L.readwith,R.readwith),
        error('lh and rh must be of the same type and format');
    end
    
    % Prepare data to save
    B = L;
    switch R.readwith,
        case {'load','csvread','fs_load_mgh'},
            B.data = cat(1,L.data,R.data);
            
        case 'dpxread',
            B.data = cat(1,L.data,R.data);
            B.extra.crd = cat(1,L.extra.crd,R.extra.crd);
            B.extra.idx = (1:size(B.data,1))';
            
        case {'srfread','fs_read_surf'},
            B.data.vtx = cat(1,L.data.vtx,R.data.vtx);
            B.data.fac = cat(1,L.data.fac,R.data.fac+size(L.data.vtx,1));
            
        case 'fs_read_curv',
            B.data = cat(1,L.data,R.data);
            B.extra.fnum = L.extra.fnum + R.extra.fnum;
            
        otherwise
            warning('Cannot deal with files read with %s. Skipping: ?%s\n',L.readwith,Flist{a}(2:end));
    end
    
    % Save
    B.filename = F{f,3};
    if any(strcmpi(B.filename(end-3:end),{'.mgh','.mgz'})),
        B.filename = B.filename(1:end-4);
    end
    palm_miscwrite(B);
end