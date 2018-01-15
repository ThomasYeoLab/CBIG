% mtimesx_build compiles mtimesx.c with BLAS libraries
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    mtimesx_build
%  Filename:    mtimesx_build.m
%  Programmer:  James Tursa
%  Version:     1.40
%  Date:        October 4, 2010
%  Copyright:   (c) 2009, 2010 by James Tursa, All Rights Reserved
%
%  This code uses the BSD License:
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
%
%--
% 
%  mtimesx_build compiles mtimesx.c and mtimesx_RealTimesReal.c with the BLAS
%  libraries libmwblas.lib (if present) or libmwlapack.lib (if libmwblas.lib
%  is not present). This function basically works as follows:
%
%  - Opens the current mexopts.bat file in the directory [prefdir], and
%    checks to make sure that the compiler selected is cl or lcc. If it
%    is not, then a warning is issued and the compilation continues with
%    the assumption that the microsoft BLAS libraries will work.
%
%  - Looks for the file libmwblas.lib or libmwlapack.lib files in the
%    appropriate directory: [matlabroot '\extern\lib\win32\microsoft']
%                       or  [matlabroot '\extern\lib\win32\lcc']
%                       or  [matlabroot '\extern\lib\win64\microsoft']
%                       or  [matlabroot '\extern\lib\win64\lcc']
%
%  - Changes directory to the directory of the file mtimesx.m.
%
%  - Compiles mtimesx.c (which includes mtimesx_RealTimesReal.c) along with
%    either libmwblas.lib or libmwlapack.lib depending on the version of
%    MATLAB. The resulting exedcutable mex file is placed in the same
%    directory as the source code. The files mtimesx.m, mtimesx.c, and
%    mtimesx_RealTimesReal.c must all be in the same directory.
%
%  - Changes the directory back to the original directory.
%
%  Change Log:
%  2009/Sep/27 --> 1.00, Initial Release
%  2010/Feb/15 --> 1.10, Fixed largearrardims typo to largeArrayDims
%  2010/Oct/04 --> 1.40, Updated support for OpenMP compiling
%
%**************************************************************************

function mtimesx_build(x)
disp(' ');
disp('... Build routine for mtimesx');

TRUE = 1;
FALSE = 0;

%\
% Check for number of inputs & outputs
%/

noopenmp = FALSE;
if( nargin == 1 )
    if( isequal(upper(x),'NOOPENMP') )
        noopenmp = TRUE;
    else
        error('Invalid input.');
    end
elseif( nargin ~= 0 )
    error('Too many inputs. Expected none.');
end
if( nargout ~= 0 )
    error('Too many outputs. Expected none.');
end

%\
% Check for non-PC
%/

disp('... Checking for PC');
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end

if( ~pc )
    disp('Non-PC auto build is not currently supported. You will have to');
    disp('manually compile the mex routine. E.g., as follows:');
    disp(' ');
    disp('>> blas_lib = ''the_actual_path_and_name_of_your_systems_BLAS_library''');
    disp('>> mex(''-DDEFINEUNIX'',''mtimesx.c'',blas_lib)');
    disp(' ');
    disp('or');
    disp(' ');
    disp('>> mex(''-DDEFINEUNIX'',''-largeArrayDims'',''mtimesx.c'',blas_lib)');
    disp(' ');
    error('Unable to compile mtimesx.c');
end

%\
% Check to see that mtimesx.c source code is present
%/

disp('... Finding path of mtimesx C source code files');
try
    mname = mfilename('fullpath');
catch
    mname = mfilename;
end
cname = [mname(1:end-6) '.c'];
if( isempty(dir(cname)) )
    disp('Cannot find the file mtimesx.c in the same directory as the');
    disp('file mtimesx_build.m. Please ensure that they are in the same');
    disp('directory and try again. The following file was not found:');
    disp(' ');
    disp(cname);
    disp(' ');
    error('Unable to compile mtimesx.c');
end
disp(['... Found file mtimesx.c in ' cname]);

%\
% Check to see that mtimesx_RealTimesReal.c source code is present
%/

rname = [mname(1:end-13) 'mtimesx_RealTimesReal.c'];
if( isempty(dir(rname)) )
    disp('Cannot find the file mtimesx_RealTimesReal.c in the same');
    disp('directory as the file mtimesx_build.m. Please ensure that');
    disp('they are in the same directory and try again. The');
    disp('following file was not found:');
    disp(' ');
    disp(rname);
    disp(' ');
    error('Unable to compile mtimesx.c');
end
disp(['... Found file mtimesx_RealTimesReal.c in ' rname]);

%\
% Open the current mexopts.bat file
%/

mexopts = [prefdir '\mexopts.bat'];
fid = fopen(mexopts);
if( fid == -1 )
    error('A C/C++ compiler has not been selected with mex -setup');
end
disp(['... Opened the mexopts.bat file in ' mexopts]);
disp('... Reading the mexopts.bat file to find the compiler and options used.');

%\
% Check for the correct compiler selected.
%/

ok_cl = FALSE;
ok_lcc = FALSE;
omp_option = '';
compiler = '(unknown)';
compilername = '';
while( TRUE )
    tline = fgets(fid);
    if( isequal(tline,-1) )
        break;
    else
        if( isempty(compilername) )
            y = findstr(tline,'OPTS.BAT');
            if( ~isempty(y) )
                x = findstr(tline,'rem ');
                if( ~isempty(x) )
                    compilername = tline(x+4:y-1);
                end
            end
        end
        x = findstr(tline,'COMPILER=lcc');
        if( ~isempty(x) )
            ok_lcc = TRUE;
            libdir = 'lcc';
            compiler = 'LCC';
            disp(['... ' compiler ' is the selected compiler']);
            break;
        end
        x = findstr(tline,'COMPILER=cl');
        if( ~isempty(x) )
            ok_cl = TRUE;
            libdir = 'microsoft';
            compiler = ['Microsoft_' compilername '_cl'];
            omp_option = ' /openmp';
            disp(['... ' compiler ' is the selected compiler']);
            break;
        end
        x = findstr(tline,'COMPILER=bcc32');
        if( ~isempty(x) )
            ok_cl = TRUE;
            libdir = 'microsoft';
            compiler = ['Borland_' compilername '_bcc32'];
            disp(['... ' compiler ' is the selected compiler']);
            disp('... Assuming that Borland will link with Microsoft libraries');
            break;
        end
        x = findstr(tline,'COMPILER=icl');
        if( ~isempty(x) )
            ok_cl = TRUE;
            if( pc )
                omp_option = ' -Qopenmp';
            else
                omp_option = ' -openmp';
            end
            libdir = 'microsoft';
            compiler = ['Intel_' compilername '_icl'];
            disp(['... ' compiler ' is the selected compiler']);
            disp('... Assuming that Intel will link with Microsoft libraries');
            break;
        end
        x = findstr(tline,'COMPILER=wc1386');
        if( ~isempty(x) )
            ok_cl = TRUE;
            libdir = 'microsoft';
            compiler = ['Watcom_' compilername '_wc1386'];
            disp(['... ' compiler ' is the selected compiler']);
            disp('... Assuming that Watcom will link with Microsoft libraries');
            break;
        end
        x = findstr(tline,'COMPILER=gcc');
        if( ~isempty(x) )
            ok_cl = TRUE;
            libdir = 'microsoft';
            omp_option = ' -fopenmp';
            compiler = 'GCC';
            disp(['... ' compiler ' is the selected compiler']);
            disp('... Assuming that GCC will link with Microsoft libraries');
            break;
        end
    end
end
fclose(fid);

%\
% MS Visual C/C++ or lcc compiler has not been selected
%/

if( ~(ok_cl | ok_lcc) )
    warning('... Supported C/C++ compiler has not been selected with mex -setup');
    warning('... Assuming that Selected Compiler will link with Microsoft libraries');
    warning('... Continuing at risk ...');
    libdir = 'microsoft';
end

%\
% If an OpenMP supported compiler is potentially present, make sure that the
% necessary compile option is present in the mexopts.bat file on the COMPFLAGS
% line.  If necessary, build a new mexopts.bat file with the correct option
% added to the COMPFLAGS line.
%/

while( TRUE )
ok_openmp = FALSE;
ok_compflags = FALSE;
xname = '';
if( isempty(omp_option) )
    disp('... OpenMP compiler not detected ... you may want to check this website:');
    disp('    http://openmp.org/wp/openmp-compilers/');
elseif( noopenmp )
    disp(['... OpenMP compiler potentially detected, but not checking for ''' omp_option ''' compile option']);
else
    disp('... OpenMP compiler potentially detected');
    disp(['... Checking to see that the ''' omp_option ''' compile option is present']);
    fid = fopen(mexopts);
    while( TRUE )
        tline = fgets(fid);
        if( isequal(tline,-1) )
            break;
        else
            x = findstr(tline,'set COMPFLAGS');
            if( ~isempty(x) )
                ok_compflags = TRUE;
                x = findstr(tline,omp_option);
                if( ~isempty(x) )
                    ok_openmp = TRUE;
                end
                break;
            end
        end
    end
    fclose(fid);
    if( ~ok_compflags )
        warning(['... COMPFLAGS line not found ... ''' omp_option ''' will not be added.']);
    elseif( ~ok_openmp )
        disp(['... The ''' omp_option ''' compile option is not present ... adding it']);
        xname = [mname(1:end-6) '_mexopts.bat'];
        disp(['... Creating custom options file ' xname ' with the ''' omp_option ''' option added.']);
        fid = fopen(mexopts);
        fidx = fopen(xname,'w');
        if( fidx == -1 )
            xname = '';
            warning(['... Unable to create custom mexopts.bat file ... ''' omp_option ''' will not be added']);
        else
            while( TRUE )
                tline = fgets(fid);
                if( isequal(tline,-1) )
                    break;
                else
                    x = findstr(tline,'set COMPFLAGS');
                    if( ~isempty(x) )
                        n = numel(tline);
                        e = n;
                        while( tline(e) < 32 )
                            e = e - 1;
                        end
                        tline = [tline(1:e) omp_option tline(e+1:n)];
                    end
                    fwrite(fidx,tline);
                end
            end
            fclose(fidx);
        end
        fclose(fid);
    end
end

%\
% Construct full file name of libmwblas.lib and libmwlapack.lib. Note that
% not all versions have both files. Earlier versions only had the lapack
% file, which contained both blas and lapack routines.
%/

comp = computer;
mext = mexext;
lc = length(comp);
lm = length(mext);
cbits = comp(max(1:lc-1):lc);
mbits = mext(max(1:lm-1):lm);
if( isequal(cbits,'64') | isequal(mbits,'64') )
    compdir = 'win64';
    largearraydims = '-largeArrayDims';
else
    compdir = 'win32';
    largearraydims = '';
end

lib_blas = [matlabroot '\extern\lib\' compdir '\' libdir '\libmwblas.lib'];
d = dir(lib_blas);
if( isempty(d) )
    disp('... BLAS library file not found, so linking with the LAPACK library');
    lib_blas = [matlabroot '\extern\lib\' compdir '\' libdir '\libmwlapack.lib'];
end
disp(['... Using BLAS library lib_blas = ''' lib_blas '''']);

%\
% Save old directory and change to source code directory
%/

cdold = cd;
if( length(mname) > 13 )
    cd(mname(1:end-13));
end

%\
% Do the compile
%/

disp('... Now attempting to compile ...');
disp(' ');
try
    if( isunix )
        if( isempty(largearraydims) )
            if( isempty(xname) )
                disp(['mex(''-DDEFINEUNIX'',''' cname ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex('-DDEFINEUNIX',cname,lib_blas,['-DCOMPILER=' compiler]);
            else
                disp(['mex(''-f'',''' xname ''',''-DDEFINEUNIX'',''' cname ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex('-f',xname,'-DDEFINEUNIX',cname,lib_blas,['-DCOMPILER=' compiler]);
            end
        else
            if( isempty(xname) )
                disp(['mex(''-DDEFINEUNIX'',''' cname ''',''' largearraydims ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex('-DDEFINEUNIX',largearraydims,cname,lib_blas,['-DCOMPILER=' compiler]);
            else
                disp(['mex(''-f'',''' xname ''',''-DDEFINEUNIX'',''' cname ''',''' largearraydims ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex('-f',xname,'-DDEFINEUNIX',largearraydims,cname,lib_blas,['-DCOMPILER=' compiler]);
            end
        end
    else
        if( isempty(largearraydims) )
            if( isempty(xname) )
                disp(['mex(''' cname ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex(cname,lib_blas,['-DCOMPILER=' compiler]);
            else
                disp(['mex(''-f'',''' xname ''',''' cname ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex('-f',xname,cname,lib_blas,['-DCOMPILER=' compiler]);
            end
        else
            if( isempty(xname) )
                disp(['mex(''' cname ''',''' largearraydims ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex(cname,largearraydims,lib_blas,['-DCOMPILER=' compiler]);
            else
                disp(['mex(''-f'',''' xname ''',''' cname ''',''' largearraydims ''',lib_blas,''-DCOMPILER=' compiler ''')']);
                disp(' ');
                mex('-f',xname,cname,largearraydims,lib_blas,['-DCOMPILER=' compiler]);
            end
        end
    end
    disp('... mex mtimesx.c build completed ... you may now use mtimesx.');
    disp(' ');
    mtimesx;
    break;
catch
    if( noopenmp )
        cd(cdold);
        disp(' ');
        disp('... Well, *that* didn''t work either!');
        disp(' ');
        disp('The mex command failed. This may be because you have already run');
        disp('mex -setup and selected a non-C compiler, such as Fortran. If this');
        disp('is the case, then rerun mex -setup and select a C/C++ compiler.');
        disp(' ');
        error('Unable to compile mtimesx.c');
    else
        disp(' ');
        disp('... Well, *that* didn''t work ...');
        disp(' ');
        if( isequal(omp_option,' /openmp') )
            disp('This may be because an OpenMP compile option was added that the');
            disp('compiler did not like. For example, the Standard versions of the');
            disp('Microsoft C/C++ compilers do not support OpenMP, only the');
            disp('Professional versions do. Attempting to compile again but this');
            disp(['time will not add the ''' omp_option ''' option.'])            
        else
            disp('This may be because an OpenMP compile option was added that the');
            disp('compiler did not like. Attempting to compile again, but this time');
            disp(['will not add the ''' omp_option ''' option.'])            
        end
        disp(' ');
        noopenmp = TRUE;
    end
end
end

%\
% Restore old directory
%/

cd(cdold);

return
end
