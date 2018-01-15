function convertmatfiles

savedir = 'matlab6files';

% take all mat files in the current directory
whdir = dir( '*.mat' );
nfiles = length( whdir );
for iiii=1:nfiles
   filenow = whdir( iiii ).name;
    
   fprintf( 'Working on file: %s\n' , filenow );
   
   % first, get the variables from thi 
   V = who( '-file' , filenow );
   
   nvars = length( V );
   vars  = [];
   for v=1:nvars
       vars = [ vars ' ' V{ v } ];
   end
   
   load( filenow );
   
   newfile = sprintf( '%s\\%s' , savedir , filenow );
   
   comm = sprintf( 'save %s %s -V6' , newfile , vars );  
   fprintf( 'comm = %s\n' , comm );
   
   eval( comm );
   
   comm = sprintf( 'clear %s' , vars );
   eval( comm );
   
end