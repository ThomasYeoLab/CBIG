
1. Overview

The nifti1 matlab i/o code is written using matlab object oriented
data structures (in matlab, see help datatypes).  The two
main structures are the nifti object, that primarily manages the
metadata (header information) and the file_array object, that manages
the data array and disk files.

The nifti1 matlab code was provided by John Ashburner, Functional
Imaging Laboratory, Wellcome Department of Imaging Neuroscience, London.  
This code is released under the GNU public license, see the license.txt 
and gpl.txt files in the distribution.
This niftimatlib release was pulled in March 2012 from the spm8  
release: "Version 4667 (SPM8) 27-Feb-12"


2. Install/Build

The nifti1 matlab i/o code was written to run under MATLAB version 6.5 or higher.
To install, just make sure that the niftimatlib/matlab directory is in your MATLAB
path.  For example. in matlab you can run the addpath command:
addpath('/usr/local/pkg/niftimatlib/matlab')
Or, you can copy the contents of the niftimatlib/matlab directory to your <home>/matlab directory.

There are two C program files included in the distribution: file2mat.c and mat2file.c
to handle file i/o.  These need to be compiled into MATLAB mex files.  
Precompiled mex files, taken from the spm8 distribution courtesy of the FIL, are included in 
this distribution for the following platforms: 
mexglx  	glnx86		Linux on x86	
mexa64  	glnxa64		Linux on x86_64
mexmaci 	maci		Apple Mac OS X on x86
mexmaci64       maci64		Apple Mac OS X on x86_64
mexw32  	win32		Microsoft Windows on x86
mexw64  	win64		Microsoft Windows on x64
So, you may not need to do the mex compile.  If you do compile,
a Makefile is in the matlab directory.  Instructions are in the Makefile, a simple
"make all" should work.  Note that you must have a MATLAB version 6.5 or higher mex compiler.
Alternately, a make.m file for calling mex from matlab was contributed by Alle Meije Wink.
Optional C code for a mex interface to Robert Cox's (NIH) nifti_stats.c code is provided
in the @nifti/private/src directory.


3. Tiny Example

Short example for those who want to see something in a hurry (longer example below):
For access to the avg152T1_LR_nifti.nii image see
http://nifti.nimh.nih.gov/nifti-1/data


To open an existing nifti1 file:

>> % be sure to rmpath any paths that point to any spm version
>> % after removing spm paths, clear all, then add niftimatlib path
>> % eg:
>> rmpath(genpath('/usr/local/pkg/spm8'))
>> clear all
>> addpath /usr/local/pkg/niftimatlib-1.2/matlab


>> f = nifti('avg152T1_LR_nifti.nii');
>> disp(f)
NIFTI object: 1-by-1
           dat: [91x109x91 file_array]
           mat: [4x4 double]
    mat_intent: 'MNI152'
          mat0: [4x4 double]
        timing: [1x1 struct]
       descrip: 'FSL3.2beta'
           cal: [0 255]
      aux_file: 'none'



>> size(f.dat)

ans =

    91   109    91


3. nifti1 i/o Class Structures


  NIFTI Object
  ------------
  
  Constructor:
  
  a = nifti(filename);
    a		- nifti object
    filename	- filename for nifti1 dataset.  Optional parameter,
    		  if omitted and empty nifti object is returned.
  

  Methods:
 
    create     - Create a NIFTI-1 file
    disp       - Disp a NIFTI-1 object
    display    - Display a NIFTI-1 object
    fieldnames - Fieldnames of a NIFTI-1 object
    nifti      - Create a NIFTI-1 object
    subsasgn   - Subscript assignment
    subsref    - Subscript referencing
    
  Fields:
      
    aux_file
    cal
    dat
    descrip
    diminfo
    intent
    mat
    mat0
    mat0_intent
    mat_intent
    timing



  FILE_ARRAY Object
  -----------------
 
  Constructor:
   
  a = file_array(fname,dim,dtype,offset,scl_slope,scl_inter)
  a         - file_array object
  fname     - filename
  dim       - dimensions (default = [0 0] )
  dtype     - datatype   (default = 'uint8-le')
  offset    - offset into file (default = 0)
  scl_slope - scalefactor (default = 1)
  scl_inter - DC offset, such that dat = raw*scale + inter (default = 0)


  Methods:

	cat		-  Concatenate file_array objects.        
	disp 		-  Display a file_array object      
	display 	-  Display a file_array object   
	double 		-  Convert to double precision     
	end 		-  Overloaded end function for file_array objects        
	fieldnames	-  Fieldnames of a file-array object 
	horzcat		-  Horizontal concatenation of file_array objects     
	length 		-  Overloaded length function for file_array objects
	ndims 		-  Number of dimensions
	numel 		-  Number of simple file arrays involved
	numeric 	-  Convert to numeric form  
	reshape 	-  Overloaded reshape function for file_array objects   
	size 		-  Overloaded size function for file_array objects     
	subsasgn 	-  Overloaded subsasgn function for file_array objects
	subsref		-  Subscripted reference.
	vertcat		-  Vertical concatenation of file_array objects
	
  Disallowed file_array methods:
	ctranspose	-  Complex conjugate transposing is not allowed
	permute		-  Permuting is not allowed
	transpose 	-  Transposing is not allowed

  Fields:

    fname
    dim
    dtype
    offset
    scl_inter
    scl_slope
	
	
	 
4. Examples

  % Example of creating a simulated .nii file.
  dat         = file_array;
  dat.fname   = 'junk.nii';
  dat.dim     = [64 64 32];
  dat.dtype   = 'FLOAT64-LE';
  dat.offset  = ceil(348/8)*8;

  % alternatively:
  % dat = file_array( 'junk.nii',dim,dtype,off,scale,inter)
  
  disp(dat)
  
  % Create an empty NIFTI structure
  N = nifti;
  
  fieldnames(N) % Dump fieldnames
  
  % Creating all the NIFTI header stuff
  N.dat = dat;
  N.mat = [2 0 0 -110 ; 0 2 0 -110; 0 0 -2 92; 0 0 0 1];
  N.mat_intent = 'xxx'; % dump possibilities
  N.mat_intent = 'Scanner';
  N.mat0 = N.mat;
  N.mat0_intent = 'Aligned';
  
  N.diminfo.slice = 3;
  N.diminfo.phase = 2;
  N.diminfo.frequency = 2;
  N.diminfo.slice_time.code='xxx'; % dump possibilities 
  N.diminfo.slice_time.code = 'sequential_increasing';
  N.diminfo.slice_time.start = 1;
  N.diminfo.slice_time.end = 32;
  N.diminfo.slice_time.duration = 3/32;
  
  N.intent.code='xxx' ; % dump possibilities
  N.intent.code='FTEST'; % or N.intent.code=4;
  N.intent.param = [4 8];
  
  N.timing.toffset = 28800;
  N.timing.tspace=3;
  N.descrip = 'This is a NIFTI-1 file';
  N.aux_file='aux-file-name.txt';
  N.cal = [0 1];
  
  create(N); % Writes hdr info
  
  %% Note that this call writes the data to the disk file
  dat(:,:,:)=0; % Write out the data as all zeros

  [i,j,k] = ndgrid(1:64,1:64,1:32);
  dat(find((i-32).^2+(j-32).^2+(k*2-32).^2 < 30^2))=1; % Write some ones in the file
  dat(find((i-32).^2+(j-32).^2+(k*2-32).^2 < 15^2))=2;
  

  % displaying a slice
  imagesc(dat(:,:,12));colorbar
  
  % get a handle to 'junk.nii';
  M=nifti('junk.nii');
  
  imagesc(M.dat(:,:,12));
 


