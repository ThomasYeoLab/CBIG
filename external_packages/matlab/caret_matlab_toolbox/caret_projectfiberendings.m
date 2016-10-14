function caret_projectfiberendings(coord_name,mat_name); 

dist_cutoff=3;

if (nargin<1 | isempty(coord_name))
    coord_name=spm_get(1,'*.coord','Pick Left Fiducial Surface');
end;
FID=caret_load(coord_name);

if (nargin<2 | isempty(mat_name))
    mat_name=spm_get(1,'*.mat','DTI - Anatomical alignment');
end;
load(mat_name);

% Make metric-file structure to save results 
Metric.encoding={'BINARY'};
Metric.num_cols=1;
Metric.num_rows=FID.num_nodes; 
Metric.column_name{1}='density';
Metric.data=zeros(Metric.num_rows,1);
Metric.column_color_mapping=[0 30];

% Load fiber-data 
load temp_results; 
num_seeds=length(fiber_stats);

dti_voxsize=[.82813 .82813 2.2];
for n=1:50 
    % Get the data into the anatomical space 
    xyzstarting_point=fiber_stats(n).xyzstarting_point./dti_voxsize+0.5; 
    xyzstarting_point=(M*[xyzstarting_point 1]')';

    xyzending_point=fiber_stats(n).xyzending_point; 
    K=size(xyzending_point,1);
    xyzending_point=xyzending_point./repmat(dti_voxsize,K,1)+0.5;
    xyzending_point=(M*[xyzending_point ones(K,1)]')';
    
    xyzending_vector=fiber_stats(n).xyzending_vector./repmat(dti_voxsize,K,1); 
    xyzending_vector=(M(1:3,1:3)*[xyzending_vector]')';
    
    % Loop over all ending fibers 
    vect_length=sqrt(sum(xyzending_vector.^2,2));
    fibers=find(vect_length>0);
    for f=1:length(fibers)
        d=FID.data-repmat(xyzending_point(fibers(f),1:3),FID.num_nodes,1);
        % make the inverse ellipsoid matrix A=
        A=eye(3);
        dist=sum(d*A.*d,2);
        [a,b]=min(dist);
        if (a<dist_cutoff)
            Metric.data(b)=Metric.data(b)+1;
        end;
        if (mod(f,10)==0)
            fprintf('.');
        end;
    end;
    fprintf('\n');
end;
keyboard;
