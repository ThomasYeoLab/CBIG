function palm_icodown(varargin)
% Downsample a data-per-vertex (DPV), data-per-face (DPF),
% surface (SRF), or surfvol (MGH/MGZ) from a higher-order tessellated
% icosahedron to a lower order one. Note that the file must have been
% derived from either FreeSurfer or Brainder tools. It may not work
% if the vertices follow a different sequence inside the file.
%
% For facewise, the downsampling method is pycnophylactic and,
% therefore, should be used for quantities that require mass
% conservation, such as areal quantities.
% For vertexwise, the method removes the redundant vertices.
% Consider smoothing if needed before downsampling.
%
% Usage:
% palm_icodown(filein,ntarget,fileout,surface)
%
% - filein  : File to be downsampled.
% - ntarget : Icosahedron order of the downsampled file.
% - fileout : Output file, downsampled.
% - suface  : Optional. For facewise, if the face indices aren't
%             ordered in the usual manner, entering a surface
%             is necessary to merge the faces correctly.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Apr/2012 (first version)
% Jul/2016 (this version, with vertexwise, surface and mgh/mgz)
% http://brainder.org

% Accept arguments
nargin  = numel(varargin);
filein  = varargin{1};
ntarget = varargin{2};
fileout = varargin{3};

% Constants for the icosahedron
V0 = 12;
F0 = 20;

% Read data from disk
X = palm_miscread(filein);
% try
%     [dpx,crd,idx] = dpxread(filein);
%     nX = numel(dpx);
%     ftype = 'dpx';
% catch
%     try
%         [vtx,fac] = srfread(filein);
%         nV = size(vtx,1);
%         ftype = 'srf';
%     catch
%         error('File %s not recognised as either DPX (curvature) or SRF (surface).',filein);
%     end
% end

switch X.readwith,
    
    case {'srfread','fs_read_surf'},
        
        % Data to operate on
        vtx = X.data.vtx;
        fac = X.data.fac;
        nV  = size(vtx,1);
        
        % Find icosahedron order
        n = round(log((nV-2)/(V0-2))/log(4));
        
        % Sanity check
        if nV ~= 4^n*(V0-2)+2,
            error('Data not from icosahedron.');
        elseif ntarget >= n,
            error('This script only downsamples data.');
        else
            fprintf('Downsampling surface geometry:');
        end
        
        % Remove vertices:
        vtx = vtx(1:(4^ntarget*(V0-2)+2),:);
        
        % Remove face indices:
        for j = (n-1):-1:ntarget,
            fprintf(' %d',j);
            nVj    = 4^j*(V0-2)+2;
            facnew = zeros(4^j*F0,3);
            fout   = find(all(fac > nVj,2));
            for f = 1:numel(fout),
                vidx = fac(fout(f),:);
                ftomerge = fac(sum(...
                    fac == vidx(1) | ...
                    fac == vidx(2) | ...
                    fac == vidx(3),2) == 2,:);
                facnew(f,:) = sum(ftomerge.*(ftomerge <= nVj),1);
            end
            fac = facnew;
        end
        fprintf('. Done.\n');
        
        % Prepare to save
        X.data.vtx = vtx;
        X.data.fac = fac;
        X.filename = fileout;
        
    case {'dpxread','fs_read_curv','fs_load_mgh'},
        
        % Data to operate on
        dpx = X.data;
        nX = size(dpx,1);
        
        % Detect what kind of data this is
        if mod(nX,10), % vertexwise
            
            % Find icosahedron order
            n = round(log((nX-2)/(V0-2))/log(4));
            
            % Sanity check
            if nX ~= 4^n*(V0-2)+2,
                error('Data not from icosahedron.');
            elseif ntarget >= n,
                error('This script only downsamples data.');
            else
                fprintf('Downsampling vertexwise data.\n');
            end
            
            % Downsample vertices
            dpx = dpx(1:(4^ntarget*(V0-2)+2),:,:,:);
            
        else % facewise
            
            % Find icosahedron order
            n = round(log(nX/F0)/log(4));
            
            % Sanity check
            if nX ~= 4^n*F0,
                error('Data not from icosahedron.');
            elseif ntarget >= n,
                error('This script only downsamples data.')
            else
                fprintf('Downsampling facewise data:');
            end
            
            % If a surface was given
            if nargin == 4,
                
                % Load it
                [~,fac] = srfread(filein);
                
                % Downsample faces (general case)
                for j = (n-1):-1:ntarget,
                    fprintf(' %d',j);
                    nVj    = 4^j*(V0-2)+2;
                    facnew = zeros(4^j*F0,3);
                    fout   = find(all(fac > nVj,2));
                    dpfnew = dpf(fout,:,:,:);
                    for f = 1:numel(fout),
                        vidx = fac(fout(f),:);
                        fidx = sum(...
                            fac == vidx(1) | ...
                            fac == vidx(2) | ...
                            fac == vidx(3),2) == 2;
                        ftomerge = fac(fidx,:);
                        facnew(f,:) = sum(ftomerge.*(ftomerge <= nVj),1);
                        dpfnew(f,:,:,:) = dpfnew(f,:,:,:) + sum(dpf(fidx,:,:,:),1);
                    end
                    fac = facnew;
                end
                fprintf('. Done.\n');
                
            else
                % Downsample faces (platonic only)
                for d = 1:(n-ntarget),
                    siz = size(dpx);
                    siz(1) = siz(1)/4;
                    dpxnew = zeros(siz);
                    for n = 1:size(dpx,4);
                        dpxnew(:,1,1,n) = sum(reshape(dpx(:,1,1,n),[4 nX/4]))';
                    end
                end
            end
        end
        
        % Prepare to save
        X.data = dpx;
        X.filename = fileout;
        switch X.readwith,
            case 'dpxread',
                nXdown = numel(dpx);
                X.extra.crd(nXdown+1:end,:) = [];
                X.extra.idx(nXdown+1:end,:) = [];
            case 'fs_read_curv',
                X.extra.fnum = 4^ntarget*F0; % this will only work for vertexwise
            case 'fs_load_mgh',
                X.extra.volsz = size(X.data);
        end
end

% Save now
palm_miscwrite(X);

