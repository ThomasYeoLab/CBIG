% svm-make.

if(ispc)
    check_dir = 'dir';
else
    check_dir = 'ls';
end

[status, y] = system([check_dir ' libsvm-mat-2.86-1']);
if(status == 0)
    cd('libsvm-mat-2.86-1');
    SVMcompile;
    cd('..');
else
   disp('libsvm-mat-2.86-1 not found, skip compilation of SVM (it is ok)'); 
end
   

[status, y] = system([check_dir ' min_heap']);
if(status == 0)
    cd('min_heap');
    min_heap_compile;
    cd('..');
else
    disp('min_heap not found, skip compilation of min_heap (it is ok)');
end

[status, y] = system([check_dir ' BasicTools']);
if(status == 0)
    cd('BasicTools');
    BasicToolsCompile;
    cd('..');
else
    disp('BasicTools not found, skip compilation of Basic Tools (it is ok)');
end
    
[status, y] = system([check_dir ' AtlasSharpnessAndLabelMatch']);
if(status == 0)
    cd('AtlasSharpnessAndLabelMatch');
    AtlasSharpnessCompile;
    cd('..');
else
    disp('AtlasSharpnessAndLabelMatch not found, skip compilation of AtlasSharpnessAndLabelMatch (it is ok)');
end
    
[status, y] = system([check_dir ' kd_tree']);
if(status == 0)
    cd('kd_tree');
    kd_tree_compile;
    cd('..');
else
    disp('kd_tree not found, skip compilation of kd_tree (it is ok)');
end

[status, y] = system([check_dir ' MARS2']);
if(status == 0)
    cd('MARS2');
    MARS2compile;
    cd('..');
else
    disp('MARS2 not found, skip compilation of MARS2 (it is ok)');
end

clear all;