@echo off

rem % SETUP command line options
rem % the first call gets us all the vc8.00 tools including lib.exe

rem % This call gives the cross-compiler for 32-bit windows to 64-bit windows
rem call D:\vstudio8\VC\bin\x86_amd64\vcvarsx86_amd64.bat

rem % This call gives the 64-bit compiler on 64-bit windows
rem call d:\vstudio8\VC\bin\amd64\vcvarsamd64.bat

echo This program must be called from the Visual Studio 2005 x64 Command Prompt
echo Make sure you have updated the path to boost

set LIBNAME=libmbgl-pcwin64-large.lib
rem OUTDIR requires the extra "\\" to avoid escaping the \" when used with quotes below
set OUTDIR=Release\\

set YASMICDIR=.
set BOOSTDIR=C:\dev\boost_1_36_0

set CFLAGS=/c /nologo /I"%YASMICDIR%" /I"%BOOSTDIR%" /Fo"%OUTDIR%" /EHsc /DNDEBUG /DMATLAB_BGL_LARGE_ARRAYS /O2 
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



