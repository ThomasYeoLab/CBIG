function rotatemat=rri_bootprocrust(origlv,bootlv)
%syntax rotatemat=rri_bootprocrust(origlv,bootlv)
% - from PLS toolbox

%define coordinate space between original and bootstrap LVs
temp=origlv'*bootlv;

%orthogonalze space
[V,W,U]=svd(temp);

%determine procrustean transform
rotatemat=U*V';
