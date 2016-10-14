% publish which files?

str.maxWidth  = 700;
str.maxHeight = 500;
str.evalCode  = false;
str.figureSnapMethod = 'print';


% publish( 'exampleLDA1' , str );
% publish( 'exampleLDA2' , str );
% publish( 'exampleVIZ1' , str );
% publish( 'exampleVIZ2' , str );
% publish( 'exampleimages1' , str );
% publish( 'exampleimages2' , str );
% publish( 'exampleAT1' , str );

str.maxWidth  = 700;
str.maxHeight = 500;
str.evalCode  = false;
str.figureSnapMethod = 'print';

%whfiles = { 'WriteTopics.m' , 'WriteTopicsMult.m' , 'VisualizeTopics.m' , 'VisualizeDocs.m' };
%whfiles = { 'OrderTopics.m' };
%whfiles = { 'importworddoccounts.m' };
whfiles = { 'GibbsSamplerAT.m' };

for j=1:length( whfiles )
    filechange = whfiles{ j };
    tempfile   = 'tempcopy.m';

    copyfile(filechange,tempfile);
    [ WL ] = textread( filechange , '%s' , 'delimiter' , '' );
    fid = fopen( filechange , 'w' );
    
    comment = true;
    for i=2:length( WL )
        if WL{i}(1) ~= '%'
            comment = false;
        end

        if comment
            fprintf( fid , '%s\r\n' , WL{ i } );
        end
    end
    
    fprintf('Publishing file %s\n' , filechange );
    publish( filechange , str );
    copyfile(tempfile,filechange);
end
