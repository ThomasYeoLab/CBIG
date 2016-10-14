opts = @(outputDir) struct(...
    'format','html',...
    'stylesheet',fullfile(pwd,'../doc/','mxdom2mbgl-html.xsl'),...
    'outputDir', fullfile(pwd,['../doc/html/' outputDir '/']));

try
    cd ('../examples');
    publish('red_black.m', opts('red_black'));
    publish('record_alg.m', opts('record_alg'));
    cnopts = opts('core_numbers_example'); cnopts.figureSnapMethod='getframe';
    publish('core_numbers_example.m', cnopts);
    publish('reweighted_graphs.m', opts('reweighted_graphs'));
    publish('new_in_3_0.m', opts('new_in_3'));
    publish('new_in_4_0.m', opts('new_in_4'));
    publish('planar_graphs.m', opts('planar_graphs'));
catch
    cd ('../doc');
end
