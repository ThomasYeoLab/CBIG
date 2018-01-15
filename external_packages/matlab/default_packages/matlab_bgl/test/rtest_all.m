D = dir;

numtests = 0;
results = [];

for di = 1:length(D)
    f = D(di);
    if f.isdir
        continue;
    end
    
    [pathstr,name,ext] = fileparts(f.name);

	if ~strcmp(ext,'.m')
		continue
    end;
    
    if strcmp('rtest_all.m', f.name)
        continue;
    end
    
    if strmatch('rtest',f.name)
        numtests = numtests+1;
        
        [pathstr,name,ext] = fileparts(f.name);
        rval = eval(name);
        results(numtests) = rval;
        
        switch(rval)
            case 1
                fprintf('%s passed\n', name);
            case 0
                fprintf('** %s FAILED!!\n', name);
        end
    end
end
