function output_args = ciftisavereset(cifti, filename, caret7command)
%Save a CIFTI file as a GIFTI external binary and then convert it to CIFTI

tic
save(cifti,[filename '.gii'],'ExternalFileBinary')
toc

%unix(['/media/1TB/matlabsharedcode/ciftiunclean.sh ' filename '.gii ' filename '_.gii']);

%unix(['mv ' filename '_.gii ' filename '.gii']);
strlength=length(filename);
if strcmp('.dscalar.nii',filename(strlength-11:strlength))
  flag=' -reset-scalars';
elseif strcmp('.dtseries.nii',filename(strlength-12:strlength))
  flag=' -reset-timepoints 1 0';
elseif strcmp('.pscalar.nii',filename(strlength-11:strlength))
  flag=' -reset-scalars';
elseif strcmp('.ptseries.nii',filename(strlength-12:strlength))
  flag=' -reset-timepoints 1 0';
else
  flag='';
end

tic
unix([caret7command ' -cifti-convert -from-gifti-ext ' filename '.gii ' filename ' ' flag]);
toc

%unix([' /bin/rm ' filename '.gii ' filename '.dat ']);
delete([filename '.gii'],[filename '.dat']);

end

