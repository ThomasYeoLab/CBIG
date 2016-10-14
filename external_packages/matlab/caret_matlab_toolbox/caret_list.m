function TabDat=caret_list(varargin);
% Function that generates a Table with clusters, local maxima 
% and corrected P-values for a certain threshold u
% function TabDat=caret_list(varargin);
% Usage: 
% Tabdat=caret_list(Surface,cSPM,u,k,varargin)
%       Surface: a surface structure (see caret_getsurface)
%       cSPM: statistical map structure (see caret_getcSPM)
%       u: height-threshold
%       k: size-threshold (mm)
%       varargin:
%           'sort_by',{'p_max'/'area'/'p_cluster'/'x_coord'}
%           'contrast',contrast-number (default is 1)
%           'sign', : Signed t-test +1 / -1 
%           'coord',avrg_PALS_file: Looks up the coordinates on a different
%                                   Surface 
% caret_list('txtlist',TabDat)
%       generates a text print-out of the table 
% 
% OUTPUT: The table has the following columns:
% c	NumN	Area	max(Z)	p(unc)	p(cor)	p(cl.)	X(mm)	Y(mm)	Z(mm)
%   c: number of cluster
%   NumN: Number of Nodes in cluster 
%   Area: area of cluster in mm
%   max(Z): maximum of statistical value in that cluster 
%   p(unc): uncorrected p-value for maximum
%   p(cor): corrected p-value for maximum
%   p(cl.): cluster-wise p-value for the cluster of that size 
%   X,Y,Z:  X-,Y-, and Z- coordinate of the maximum in the space the map is
%           in. Only if the map is based on a fiducial in a MNI-space do these 
%           correspond to MNI-coordinates  
%-----------------------------------------------------------------------
% v.1.0 Joern Diedrichsen 12/10/04 jdiedric@bme.jhu.edu
if (~ischar(varargin{1}))
    S=varargin{1};
    cSPM=varargin{2};
    u=varargin{3};
    k=varargin{4};
    c=5; 
    sort_by='p_max';
    contrast=1;
    sign=1;
    coord=[];
    while c<=length(varargin)
        switch (varargin{c})
            case {'sort_by','contrast','sign','coord'}
                eval([varargin{c} '=varargin{c+1};']);
                c=c+2;
            otherwise
                error('unknown option');
        end;
    end;
    if ~isfield(S,'A');
        fprintf('Calculate Search Area\n');
        S=caret_calcarea(S);
    end;
    if ~isfield(cSPM,'FWHM');
        fprintf('Calculate Smoothness\n');
        cSPM.FWHM=caret_estsmoothness(cSPM,S);
    end;
    if ~isfield(S.Nodes,'Neighbor')
        S=caret_calcneighbor(S);
    end;
    D=length(cSPM.FWHM);
    fwhm=(prod(cSPM.FWHM).^(1./D));       % We will do stuff in mm
    x=[0:2];
    R=S.A./((fwhm).^x);
    k=k./(fwhm^2);    % K is now also in resel  
    
    if (sign==-1)
        cSPM.con(contrast).Z=-cSPM.con(contrast).Z;
    end;

    if (~isempty(coord))
        C=caret_load(coord);
    end;
    
    %-----------------------------------------------------------------------
    % Calculate Clusters and local maxima, restrict to subset if surface
    % subset is give 
    con=cSPM.con(contrast);
    if (isfield(S,'sub_index'))
        con.Z=con.Z(S.sub_index,:);    
    end;
    A=caret_clusters(S,con.Z>u);
    a=unique(A);
    cluster=0;
    TabDat.dat=[];
    for j=1:length(a)-1
        i=find(A==j);
        area=sum(S.Nodes.area(i));
        if (area./(fwhm^2))>k 
            cluster=cluster+1;
            TabDat.dat(cluster,1)=cluster;
            TabDat.dat(cluster,2)=length(i);
            TabDat.dat(cluster,3)=area;
            [maxZ,maxInd]=max(con.Z(i));
            TabDat.dat(cluster,4)=maxZ;
            TabDat.dat(cluster,5)=caret_P(1,0,maxZ,con.df,con.STAT,1);
            TabDat.dat(cluster,6)=caret_P(1,0,maxZ,con.df,con.STAT,R);
            TabDat.dat(cluster,7)=caret_P(1,TabDat.dat(cluster,3)/(fwhm^2),u,con.df,con.STAT,R);
            if (~isempty(coord))
                if (isfield(S,'sub_index'))
                    TabDat.dat(cluster,8:10)=C.data(S.sub_index(i(maxInd)),1:3);
                else
                    TabDat.dat(cluster,8:10)=C.data(i(maxInd),:);
                end;
            else
                TabDat.dat(cluster,8:10)=S.Nodes.data(i(maxInd),:);
            end;
        end;
    end;
    Pz              = caret_P(1,0,u,con.df,con.STAT,1);
    Pu              = caret_P(1,0,u,con.df,con.STAT,R);
    [P Pn Em En EN] = caret_P(1,k,u,con.df,con.STAT,R);
    
    %[P,p,Em,En,EN] = caret_P(num_clusters,4,u,cSPM.df,cSPM.STAT,R);
    %-----------------------------------------------------------------------
    %-Headers for text table...
    TabDat.tit = cSPM.title;
    TabDat.hdr = {'  c','NumN','Area','max(Z)','p(unc)','p(cor)','p(cl.)','X(mm)','Y(mm)','Z(mm)'};
    TabDat.fmt = {'%3d','%4d','%-3.2f','%0.3f','%0.3f', '%0.3f', '%0.3f','%4.2f','%4.2f','%4.2f'};				%-XYZ
    
    %-Footnote with SPM parameters
    %-----------------------------------------------------------------------
    TabDat.ftr    = cell(4,2);
    TabDat.ftr{1} = ...
        sprintf('Height threshold: %c = %0.2f, p = %0.3f (%0.3f)',...
        con.STAT,u,Pz,Pu);
    TabDat.ftr{2} = ...
        sprintf('Extent threshold: k = %0.2f sqmm, p = %0.3f (%0.3f)',...
        k*fwhm^2,Pn,P);
    TabDat.ftr{3} = ...
        sprintf('Expected sqmm per cluster, <k> = %0.3f',En*fwhm.^2);
    TabDat.ftr{4} = ...
        sprintf('Expected number of clusters, <c> = %0.2f',Em*Pn);
    TabDat.ftr{5} = ...
        sprintf('Expected surface above threshold (sqmm) = %0.2f',Em*En*fwhm.^2);
    TabDat.ftr{6} = ...
        sprintf('Degrees of freedom = [%0.1f, %0.1f]',con.df);
    TabDat.ftr{7} = ...
        sprintf('Smoothness FWHM = %0.1f (mm) ',fwhm);
    TabDat.ftr{8} = ...
        sprintf('Search vol: %5.0f sqmm; %0.1f resels',S.A(end),R(end));
    
    % -----------------------------------------------------------------
    % Put all the necessary information in Table-strcuture 
    TabDat.Em=Em;
    TabDat.En=En*fwhm.^2;
    TabDat.Pn=Pn;
    TabDat.EN=EN*fwhm.^2;
    if (isempty(TabDat.dat))
    TabDat.TotA=0;
    TabDat.m=0;
    TabDat.n=0;
    else
    TabDat.TotA=sum(TabDat.dat(:,3));
    TabDat.m=cluster;
    TabDat.n=mean(TabDat.dat(:,3));
    % -----------------------------------------------------------------
    % Sort and list the Table  
    switch (sort_by)
        case 'p_max'
            sort_col=-TabDat.dat(:,4);
        case 'area'
            sort_col=-TabDat.dat(:,3);
        case 'p_cluster'
            sort_col=TabDat.dat(:,7);
        case 'x_coord'
            sort_col=TabDat.dat(:,8);
    end;
    [A,IND]=sort(sort_col);
    TabDat.dat=TabDat.dat(IND,:);
    caret_list('txtlist',TabDat);
    end;
else 
    switch(varargin{1})
        % -----------------------------------------------------------------
        % -Print ASCII text table
        case 'txtlist'           
            TabDat = varargin{2};
            
            %-Table Title
            %-----------------------------------------------------------------------
            fprintf('\n\nSTATISTICS: %s\n',TabDat.tit);
            fprintf('%c','='*ones(1,80)), fprintf('\n');
            
            %-Table header
            %-----------------------------------------------------------------------
            fprintf('%s\t',TabDat.hdr{1:end-1});fprintf('%s\n',TabDat.hdr{end});
            fprintf('%c','-'*ones(1,80)), fprintf('\n');
            
            
            %-Table data
            %-----------------------------------------------------------------------
            for i = 1:size(TabDat.dat,1)
                for j=1:size(TabDat.dat,2)
                    fprintf(TabDat.fmt{j},TabDat.dat(i,j))
                    fprintf('\t')
                end
                fprintf('\n')
            end
            %      fprintf('%s\n',TabDat.str)
            fprintf('%c','-'*ones(1,80)), fprintf('\n')
            
            %-Table footer
            %-----------------------------------------------------------------------
            fprintf('%s\n',TabDat.ftr{:})
            fprintf('%c','='*ones(1,80)), fprintf('\n\n')
    end;
end;