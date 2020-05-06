## Script copies basic freesurfer files
##Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

for i in fsaverage fsaverage5 fsaverage6;
do mkdir ../FreeSurfer/$i/label/
mkdir ../FreeSurfer/$i/surf/
rsync -az $FREESURFER_HOME/subjects/${i}/label/*cortex.label ../FreeSurfer/$i/label/
rsync -az $FREESURFER_HOME/subjects/${i}/label/*Medial*.label ../FreeSurfer/$i/label/
rsync -az $FREESURFER_HOME/subjects/${i}/label/*aparc* ../FreeSurfer/$i/label/
for j in white orig pial inflated curv sulc;
do
rsync -az $FREESURFER_HOME/subjects/${i}/surf/*${j} ../FreeSurfer/$i/surf/
done;
done;

