@echo off

rem % SETUP command line options
rem % the first call gets us all the vc8.00 tools including lib.exe
rem % the second call gets us all the vc7.00 tools for linking with matlab

call d:\vstudio8\Common7\Tools\vsvars32.bat
rem call d:\vstudio8\vctoolkit\vcvars32.bat

set LIBNAME=libmbgl-pcwin32.lib
set OUTDIR=Release\\

set YASMICDIR=.
set BOOSTDIR=e:\dev\lib\boost_1_36_0

set CFLAGS=-c -nologo -I"%YASMICDIR%" -I"%BOOSTDIR%" /Fo"%OUTDIR%\\" /EHsc /D "NDEBUG" /O2 /ML
rem set CFLAGS=-c -nologo -Ie:\dev\yasmic -Ie:\dev\lib\boost_1_33_1 /Fo"%OUTDIR%\\" /EHsc /ML /Od /D "_DEBUG" /Fd"%OUTDIR%\vc70.pdb" /Zi
set LIBFLAGS=-nologo /out:"%OUTDIR%\\%LIBNAME%"

rem Make sure the release directory exists
echo "Checking for the release directory..."
if not exist Release ( mkdir Release ) else ( echo "Release direcotry exists" ) 


@echo on

cl %CFLAGS% components.cc
cl %CFLAGS% max_flow.cc
cl %CFLAGS% orderings.cc
cl %CFLAGS% searches.cc
cl %CFLAGS% shortest_path.cc
cl %CFLAGS% spanning_trees.cc
cl %CFLAGS% statistics.cc
cl %CFLAGS% layouts.cc
cl %CFLAGS% planar.cc

lib %LIBFLAGS% ^
  %OUTDIR%\components.obj ^
  %OUTDIR%\max_flow.obj ^
  %OUTDIR%\orderings.obj ^
  %OUTDIR%\searches.obj ^
  %OUTDIR%\shortest_path.obj ^
  %OUTDIR%\spanning_trees.obj ^
  %OUTDIR%\statistics.obj ^
  %OUTDIR%\layouts.obj ^
  %OUTDIR%\planar.obj 


